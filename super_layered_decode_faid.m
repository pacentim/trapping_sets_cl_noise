function [hard_bits, success, iters, Lpost_out, state] = super_layered_decode_faid(G, layers, syndrome, priors, cfg, state)
%SUPER_LAYERED_DECODE_FAID Layered FAID decoder with 3-bit (7-level) messages.
%
% Alphabet: {-3,-2,-1,0,1,2,3} (stored as int8)
%
% CN update: quantized min-sum (layered)
% VN update: FAID LUT for dv=2 and dv=3 (fallback to quantized sum for other dv)
%
% You will implement LUT2/LUT3 later (placeholders included).



    % -------- defaults --------
    if nargin < 5 || isempty(cfg), cfg = struct(); end
    cfg = apply_superlayered_defaults(cfg);

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

    % =======================
    % FAID tables (PLACEHOLDER)
    % =======================
    % Convention:
    %   levels = [-3 -2 -1 0 1 2 3]
    %   indices use off = 4 so that idx = val + off maps:
    %      -3->1, -2->2, -1->3, 0->4, 1->5, 2->6, 3->7
    %
    % LUT2: out = f(y, m)            where dv=2, m is the OTHER incoming c2v (excluding target edge)
    % LUT3: out = f(y, m1, m2)       where dv=3, (m1,m2) are the two OTHER incoming c2v
    %
    % Fill these later with your optimized FAID rule.
    levels = int8(-3:3);
    off    = int8(4);

    % dv = 2 (one incoming message): vector
    % LUT2neg(i) = Φv(-C, Mi)   where Mi is the incoming message level
    levels = int8(-3:3);


    LUT3 = int8([ ...
        -3 -3 -3 -3 -3 -3 -1;  % m1=-L3:  -L3 -L3 -L3 -L3 -L3 -L3 -L1
        -3 -3 -3 -3 -2 -1  1;  % m1=-L2:  -L3 -L3 -L3 -L3 -L2 -L1 +L1
        -3 -3 -2 -2 -1 -1  1;  % m1=-L1:  -L3 -L3 -L2 -L2 -L1 -L1 +L1
        -3 -3 -2 -1  0  0  1;  % m1=0:    -L3 -L3 -L2 -L1  0   0  +L1
        -3 -2 -1  0  0  1  2;  % m1=+L1:  -L3 -L2 -L1  0   0  +L1 +L2
        -3 -1 -1  0  1  1  3;  % m1=+L2:  -L3 -L1 -L1  0  +L1 +L1 +L3
        -1  1  1  1  2  3  3   % m1=+L3:  -L1 +L1 +L1 +L1 +L2 +L3 +L3
        ]);


    % -------- quantize priors --------
    % If you later want scaling, set cfg.prior_scale (e.g., 1/step) and do:
    % yq = q7(round(priors(v)*cfg.prior_scale));
    if ~isfield(cfg,'prior_scale'), cfg.prior_scale = 0.2; end
    yq = q7(round(priors * cfg.prior_scale));     % int8(N,1)

    % -------- message initialization --------
    % store messages as int8 in {-3..3}
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

                % min1/min2 and argmin
                [min1, min_idx] = min(abs_msgs);
                if numel(abs_msgs) >= 2
                    abs_msgs2 = abs_msgs;
                    abs_msgs2(min_idx) = intmax('int16');
                    min2 = min(abs_msgs2);
                else
                    min2 = intmax('int16'); % degree-1 CN case
                end

                % sign product (with signnz: treat 0 as +)
                sgn = ones(size(msgs),'int8');
                sgn(msgs < 0) = int8(-1);
                sign_prod = int8(prod(int16(sgn))); % safe multiply in int16 then cast

                if syndrome(c) ~= 0
                    sign_prod = -sign_prod;
                end

                % outgoing messages
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

                    % beta_code support (quantized)
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

                yv = yq(v); % int8, represents channel sign (+ or -)

                if dv == 2
                    e1 = edges_v(1); e2 = edges_v(2);
                    m1 = int8(G.c2v(e1));
                    m2 = int8(G.c2v(e2));

                    out1 = faid_vn_dv2(yv, m2); % msg on e1 uses other incoming m2
                    out2 = faid_vn_dv2(yv, m1); % msg on e2 uses other incoming m1

                    % optional damping (quantize after)
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

                        out = faid_vn_dv3(yv, m1, m2, LUT3);

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
            Lpost_out(v) = double(Lpost_q);          % still returning as double for compatibility
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
    if ~isfield(cfg,'max_iter'),       cfg.max_iter = 50; end
    if ~isfield(cfg,'beta_code'),      cfg.beta_code = 1.0; end
    if ~isfield(cfg,'Tmin'),           cfg.Tmin = 1.0; end
    if ~isfield(cfg,'layer_damping'),  cfg.layer_damping = 1.0; end  % for FAID you might prefer 1.0 initially
    if ~isfield(cfg,'enable_gumbel'),  cfg.enable_gumbel = false; end
end

function xq = q7(x)
% Quantize to 7 levels {-3,-2,-1,0,1,2,3}, output int8.
% Accepts double/int16/int32; rounds if not integer.
    if ~isinteger(x)
        x = round(x);
    end
    x = max(min(x, 3), -3);
    xq = int8(x);
end

function touched = precompute_layer_vars(G, layers)
% touched{l} contains unique VNs appearing in layer checks.
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
% Same as your previous Gumbel-based layer shuffling.
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


function out = faid_vn_dv2(y, m)
% FAID_VN_DV2 for degree-2
% This implementation respects the specific 'C' (magnitude) of the bit 'y'.
% If y is L1 (+1), it adds 1. If y is L3 (+3), it adds 3.
% If y is L0 (0), it adds 0.
    
    % Direct Linear Threshold: Q(m + y)
    val = int16(y) + int16(m);
    out = q7(val); 
end


function out = faid_vn_dv3(y, m1, m2, LUT3neg)
% y, m1, m2 are int8 in {-3..3}
% LUT3neg is the Table-II LUT for -C

    off = int16(4); % maps -3..3 -> 1..7

    if y < 0
        out = LUT3neg(int16(m1)+off, int16(m2)+off);
    else
        out = -LUT3neg(int16(-m1)+off, int16(-m2)+off);
    end
end

