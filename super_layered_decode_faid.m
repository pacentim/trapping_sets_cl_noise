function [hard_bits, success, iters, Lpost_out, state] = super_layered_decode_faid(G, layers, syndrome, priors, cfg, state, lut)
%SUPER_LAYERED_DECODE_FAID Layered FAID decoder with 3-bit (7-level) messages.
%
% Alphabet: {-3,-2,-1,0,1,2,3} (stored as int8)
%
% CN update: quantized min-sum (layered)
% VN update: FAID LUT for dv=2 and dv=3
%
% Optional input "lut":
%   - If omitted or empty: uses built-in default LUTs.
%   - If provided: struct with fields (any subset allowed)
%       lut.LUT2neg : int8(1x7)  dv=2 canonical table for y<0, indexed by m=-3..3
%       lut.LUT3neg : int8(7x7)  dv=3 canonical table for y<0, indexed by m1,m2=-3..3
%
% Symmetry convention (as in your dv=3):
%   if y<0:  out = LUTneg(m,...) 
%   if y>0:  out = -LUTneg(-m,...)

    % -------- defaults --------
    if nargin < 5 || isempty(cfg), cfg = struct(); end
    cfg = apply_superlayered_defaults(cfg);

    % -------- lut defaults --------
    if nargin < 7 || isempty(lut)
        lut = struct();
    end
    if ~isfield(lut,'LUT2neg') || isempty(lut.LUT2neg)
        % Default dv=2: linear rule Q(y+m) implemented as a LUT for canonical y<0.
        % For dv=2, table is indexed only by m; y is handled by symmetry and/or direct addition outside the LUT.
        % Here we set it to the identity (out = m) for y<0 canonical, and symmetry will handle y>0.
        % If you prefer a different "standard", replace this vector.
        lut.LUT2neg = int8([-3 -2 -1 0 1 2 3]);
    else
        lut.LUT2neg = int8(lut.LUT2neg(:)).'; % force row
    end
    if numel(lut.LUT2neg) ~= 7
        error('lut.LUT2neg must be 1x7 (or 7x1) int8 for m=-3..3.');
    end

    if ~isfield(lut,'LUT3neg') || isempty(lut.LUT3neg)
        % Default dv=3 LUT (your current one), canonical for y<0.
        lut.LUT3neg = int8([ ...
            -3 -3 -3 -3 -3 -3 -1;
            -3 -3 -3 -3 -2 -1  1;
            -3 -3 -2 -2 -1 -1  1;
            -3 -3 -2 -1  0  0  1;
            -3 -2 -1  0  0  1  2;
            -3 -1 -1  0  1  1  3;
            -1  1  1  1  2  3  3
        ]);
       % lut.LUT3neg = int8([ -3 -3 -3 -3 -3 -2 -1; -3 -3 -3 -3 -2 -1 0; -3 -3 -3 -2 -1 0 1; -3 -3 -2 -1 0 1 2; -3 -2 -1 0 1 2 3; -2 -1 0 1 2 3 3;  -1 0 1 2 3 3 3 ]);
    else
        lut.LUT3neg = int8(lut.LUT3neg);
    end
    if ~isequal(size(lut.LUT3neg), [7 7])
        error('lut.LUT3neg must be 7x7 int8 indexed by m1,m2=-3..3.');
    end

    % -------- basic checks --------
    M = G.M; N = G.N;
    syndrome = syndrome(:);
    priors   = priors(:);

    if numel(syndrome) ~= M
        error('Syndrome size mismatch: must match graph rows (M).');
    end
    if numel(priors) ~= N
        error('Priors size mismatch: must equal graph.N.');
    end

    L = numel(layers);

    % -------- init / precompute state --------
    if nargin < 6 || isempty(state)
        state = struct();
    end

    if ~isfield(state, 'layer_touched_vars') || isempty(state.layer_touched_vars)
        state.layer_touched_vars = precompute_layer_vars(G, layers);
    end
    if ~isfield(state, 'layer_order') || isempty(state.layer_order)
        state.layer_order = 1:L;
    else
        state.layer_order = 1:L; % mimic C++: reset each call
    end
    if ~isfield(state, 'probabilities') || numel(state.probabilities) ~= L
        state.probabilities = zeros(L,1);
    end
    if ~isfield(state, 'score_buffer') || size(state.score_buffer,1) ~= L
        state.score_buffer = zeros(L,2);
    end
    if ~isfield(state, 'rng') || isempty(state.rng)
        state.rng = RandStream('mt19937ar','Seed','shuffle');
    end

    layer_max_residual = zeros(L,1);

    levels = int8(-3:3);
    off    = int16(4); %#ok<NASGU> % maps -3..3 -> 1..7

    % -------- quantize priors --------
    if ~isfield(cfg,'prior_scale'), cfg.prior_scale = 0.5; end
    yq = q7(round(priors * cfg.prior_scale));     % int8(N,1)

    % -------- message initialization --------
    for v = 1:N
        edges_v = G.vn_edge{v};
        G.v2c(edges_v) = yq(v);
        G.c2v(edges_v) = int8(0);
    end

    hard_bits = zeros(N,1,'uint8');
    Lpost_out = zeros(N,1);
    success = false;

    for iter = 1:cfg.max_iter
        layer_max_residual(:) = 0.0;

        % ---- layered sweep ----
        for ord_pos = 1:L
            l_idx = state.layer_order(ord_pos);
            current_checks = layers{l_idx};

            max_res = 0.0;

            % ======================
            % CN update (quantized min-sum)
            % ======================
            for cc = 1:numel(current_checks)
                c = current_checks(cc);
                edges = G.cn_edge{c};
                if isempty(edges), continue; end

                msgs = int8(G.v2c(edges));               % int8
                abs_msgs = abs(int16(msgs));             % int16 for safe abs

                [min1, min_idx] = min(abs_msgs);
                if numel(abs_msgs) >= 2
                    abs_msgs2 = abs_msgs;
                    abs_msgs2(min_idx) = intmax('int16');
                    min2 = min(abs_msgs2);
                else
                    min2 = intmax('int16');
                end

                sgn = ones(size(msgs),'int8');
                sgn(msgs < 0) = int8(-1);
                sign_prod = int8(prod(int16(sgn)));

                if syndrome(c) ~= 0
                    sign_prod = -sign_prod;
                end

                for k = 1:numel(edges)
                    if k == min_idx
                        mag = min2;
                    else
                        mag = min1;
                    end

                    self_sign = int8(1);
                    if G.v2c(edges(k)) < 0
                        self_sign = int8(-1);
                    end
                    out_sign = int16(sign_prod) * int16(self_sign);

                    new_val = double(out_sign) * double(mag) * cfg.beta_code;
                    new_msg = q7(round(new_val));

                    old_msg = int8(G.c2v(edges(k)));
                    diff = abs(double(new_msg) - double(old_msg));
                    if diff > max_res
                        max_res = diff;
                    end

                    G.c2v(edges(k)) = new_msg;
                end
            end

            layer_max_residual(l_idx) = max_res;

            % ======================
            % VN updates (FAID LUT for dv=2,3)
            % ======================
            touched = state.layer_touched_vars{l_idx};
            for tt = 1:numel(touched)
                v = touched(tt);
                edges_v = G.vn_edge{v};
                dv = numel(edges_v);

                yv = yq(v); % int8 in {-3..3}

                if dv == 2
                    e1 = edges_v(1); e2 = edges_v(2);
                    m1 = int8(G.c2v(e1));
                    m2 = int8(G.c2v(e2));

                    out1 = faid_vn_dv2_lut(yv, m2, lut.LUT2neg); % msg on e1 uses other incoming m2
                    out2 = faid_vn_dv2_lut(yv, m1, lut.LUT2neg); % msg on e2 uses other incoming m1

                    if cfg.layer_damping < 1.0
                        out1 = q7(round(cfg.layer_damping*double(out1) + (1.0-cfg.layer_damping)*double(G.v2c(e1))));
                        out2 = q7(round(cfg.layer_damping*double(out2) + (1.0-cfg.layer_damping)*double(G.v2c(e2))));
                    end

                    G.v2c(e1) = out1;
                    G.v2c(e2) = out2;

                elseif dv == 3
                    e = edges_v(:);
                    m = int8(G.c2v(e)); % 3x1

                    for ii = 1:3
                        idx_other = setdiff(1:3, ii);
                        m1 = m(idx_other(1));
                        m2 = m(idx_other(2));

                        out = faid_vn_dv3_lut(yv, m1, m2, lut.LUT3neg);

                        if cfg.layer_damping < 1.0
                            out = q7(round(cfg.layer_damping*double(out) + (1.0-cfg.layer_damping)*double(G.v2c(e(ii)))));
                        end

                        G.v2c(e(ii)) = out;
                    end
                end
            end
        end

        % ---- hard decision + Lpost ----
        for v = 1:N
            edges_v = G.vn_edge{v};
            Lpost_q = int16(yq(v)) + sum(int16(G.c2v(edges_v)));
            Lpost_out(v) = double(Lpost_q);
            hard_bits(v) = uint8(Lpost_q < 0);
        end

        % ---- syndrome check ----
        success = true;
        for c = 1:M
            vs = G.cn_to_vn{c};
            parity = uint8(0);
            for kk = 1:numel(vs)
                parity = bitxor(parity, hard_bits(vs(kk)));
            end
            if parity ~= uint8(syndrome(c))
                success = false;
                break;
            end
        end

        if success
            iters = iter;
            state.G = G;
            return;
        end

        if cfg.enable_gumbel
            state.layer_order = update_layer_order(layer_max_residual, cfg, iter, state.rng);
        end
    end

    iters = cfg.max_iter;
    state.G = G;
end

% ================= helpers =================

function cfg = apply_superlayered_defaults(cfg)
    if ~isfield(cfg,'max_iter'),       cfg.max_iter = 10; end
    if ~isfield(cfg,'beta_code'),      cfg.beta_code = 1.0; end
    if ~isfield(cfg,'Tmin'),           cfg.Tmin = 1.0; end
    if ~isfield(cfg,'layer_damping'),  cfg.layer_damping = 1.0; end
    if ~isfield(cfg,'enable_gumbel'),  cfg.enable_gumbel = true; end
end

function xq = q7(x)
    if ~isinteger(x)
        x = round(x);
    end
    x = max(min(x, 3), -3);
    xq = int8(x);
end

function touched = precompute_layer_vars(G, layers)
    N = G.N;
    L = numel(layers);
    touched = cell(L,1);
    seen = false(N,1);

    for l = 1:L
        vars = [];
        checks = layers{l};
        for cc = 1:numel(checks)
            c = checks(cc);
            vs = G.cn_to_vn{c};
            for kk = 1:numel(vs)
                v = vs(kk);
                if ~seen(v)
                    seen(v) = true;
                    vars(end+1) = v; %#ok<AGROW>
                end
            end
        end
        touched{l} = vars;
        seen(vars) = false;
    end
end

function order = update_layer_order(layer_max_residual, cfg, k, rng)
    L = numel(layer_max_residual);

    pos = layer_max_residual(layer_max_residual > 0);
    if isempty(pos)
        mean_pos = 0.0;
    else
        mean_pos = mean(pos);
    end

    cooling = 1.0 - double(k - 1) / (1.2 * cfg.max_iter);
    Tk = max(mean_pos * cooling, cfg.Tmin);
    Tk = max(Tk, cfg.Tmin);

    z = layer_max_residual / Tk;
    max_z = max(z);
    probs = exp(z - max_z);
    sum_exp = sum(probs);

    if sum_exp == 0 || isinf(sum_exp)
        probs = ones(L,1) / L;
    else
        probs = probs / sum_exp;
    end

    u = rand(rng, L, 1);
    u = max(u, eps);
    u = min(u, 1 - eps);
    g = -log(-log(u));

    score = log(probs + 1e-20) + g;
    [~, perm] = sort(score, 'descend');
    order = perm(:).';
end

function out = faid_vn_dv2_lut(y, m, LUT2neg)
% LUT2neg: 1x7 canonical for y<0, indexed by m=-3..3.
    off = int16(4);
    if y < 0
        out = LUT2neg(int16(m) + off);
    else
        out = -LUT2neg(int16(-m) + off);
    end
end

function out = faid_vn_dv3_lut(y, m1, m2, LUT3neg)
% LUT3neg: 7x7 canonical for y<0, indexed by m1,m2=-3..3.
    off = int16(4);
    if y < 0
        out = LUT3neg(int16(m1)+off, int16(m2)+off);
    else
        out = -LUT3neg(int16(-m1)+off, int16(-m2)+off);
        % (your convention: y>0 -> negate and index at (-m1,-m2))
    end
end
