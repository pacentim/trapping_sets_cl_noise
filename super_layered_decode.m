function [hard_bits, success, iters, Lpost_out, state] = super_layered_decode(G, layers, syndrome, priors, cfg, state)
%SUPER_LAYERED_DECODE Layered min-sum decoder (single graph), MATLAB port.
%
% Inputs:
%   G        : Tanner graph struct from build_tanner_graph
%   layers   : cell(L,1) or {L} where layers{l} is a vector of CN indices (1..M)
%   syndrome : (M x 1) or (1 x M) vector in {0,1}
%   priors   : (N x 1) or (1 x N) double, LLRs
%   cfg      : struct with fields:
%              max_iter, saturation, beta_code, Tmin, layer_damping, enable_gumbel
%   state    : optional persistent state struct:
%              layer_touched_vars, rng, layer_order, probabilities, score_buffer
%
% Outputs:
%   hard_bits : (N x 1) uint8 (0/1)
%   success   : logical
%   iters     : number of iterations executed
%   Lpost_out : (N x 1) double
%   state     : updated state (for reuse across calls)

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
        state.score_buffer = zeros(L,2); % [score, idx]
    end
    if ~isfield(state, 'rng') || isempty(state.rng)
        state.rng = RandStream('mt19937ar','Seed','shuffle');
    end

    layer_max_residual = zeros(L,1);

    % -------- message initialization (same as C++) --------
    sat = cfg.saturation;

    for v = 1:N
        val = saturate(priors(v), sat);
        edges_v = G.vn_edge{v};
        G.v2c(edges_v) = val;
        G.c2v(edges_v) = 0.0;
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

            % CN update (min-sum with syndrome)
            for cc = 1:numel(current_checks)
                c = current_checks(cc);
                edges = G.cn_edge{c};
                if isempty(edges), continue; end

                msgs = G.v2c(edges);
                abs_msgs = abs(msgs);

                % min1/min2 and argmin
                [min1, min_idx] = min(abs_msgs);
                if numel(abs_msgs) >= 2
                    abs_msgs2 = abs_msgs;
                    abs_msgs2(min_idx) = inf;
                    min2 = min(abs_msgs2);
                else
                    min2 = inf; % degree-1 CN case
                end

                % sign product (with signnz)
                sgn = ones(size(msgs));
                sgn(msgs < 0) = -1;
                sign_prod = prod(sgn);

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

                    self_sign = 1;
                    if G.v2c(edges(k)) < 0
                        self_sign = -1;
                    end
                    out_sign = sign_prod * self_sign;

                    new_msg = saturate(out_sign * mag * cfg.beta_code, sat);
                    old_msg = G.c2v(edges(k));

                    diff = abs(new_msg - old_msg);
                    if diff > max_res
                        max_res = diff;
                    end

                    G.c2v(edges(k)) = new_msg;
                end
            end

            layer_max_residual(l_idx) = max_res;

            % VN updates on touched vars
            touched = state.layer_touched_vars{l_idx};
            for tt = 1:numel(touched)
                v = touched(tt);
                edges_v = G.vn_edge{v};
                sum_c2v = sum(G.c2v(edges_v));
                Lpost = priors(v) + sum_c2v;

                for ee = 1:numel(edges_v)
                    e = edges_v(ee);
                    extrinsic = Lpost - G.c2v(e);
                    prev = G.v2c(e);
                    out = cfg.layer_damping * extrinsic + (1.0 - cfg.layer_damping) * prev;
                    G.v2c(e) = saturate(out, sat);
                end
            end
        end

        % ---- hard decision + Lpost ----
        for v = 1:N
            edges_v = G.vn_edge{v};
            Lpost = priors(v) + sum(G.c2v(edges_v));
            Lpost_out(v) = Lpost;
            hard_bits(v) = uint8(Lpost < 0);
        end

        % ---- syndrome check ----
        success = true;
        for c = 1:M
            vs = G.cn_to_vn{c};
            parity = 0;
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
            state.G = G; % optionally return updated messages
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
    if ~isfield(cfg,'saturation'),     cfg.saturation = 25.0; end
    if ~isfield(cfg,'beta_code'),      cfg.beta_code = 1.0; end
    if ~isfield(cfg,'Tmin'),           cfg.Tmin = 1.0; end
    if ~isfield(cfg,'layer_damping'),  cfg.layer_damping = 0.9; end
    if ~isfield(cfg,'enable_gumbel'),  cfg.enable_gumbel = true; end
end

function y = saturate(x, limit)
    y = min(max(x, -limit), limit);
end

function touched = precompute_layer_vars(G, layers)
% touched{l} contains unique VNs appearing in layer checks, like C++.
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
% Port of updateLayerOrder(int k) in C++.
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

    % Gumbel sampling + sorting by score
    u = rand(rng, L, 1);
    u = max(u, eps);
    u = min(u, 1 - eps);
    g = -log(-log(u)); % gumbel

    score = log(probs + 1e-20) + g;

    [~, perm] = sort(score, 'descend');
    order = perm(:).'; % row vector of layer indices
end
