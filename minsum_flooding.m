function [hard_bits, success, iters, Lpost_out, G] = minsum_flooding(G, syndrome, priors, cfg)
% Flooding normalized min-sum on the full graph (no degree-1 injection rule).
%
% Inputs:
%   G: graph struct with fields
%       - M, N
%       - vn_edge{v} : edge indices adjacent to VN v
%       - cn_edge{c} : edge indices adjacent to CN c
%       - v2c(edge), c2v(edge) : message arrays (double)
%       - cn_to_vn{c} : VN indices adjacent to CN c (for syndrome check)
%   syndrome (M x 1) : 0/1
%   priors   (N x 1) : channel LLRs
%   cfg: optional struct with fields
%       - max_iter (default 50)
%       - saturation (default 25.0)
%       - beta_code (default 1.0)  (normalized min-sum factor)
%
% Outputs:
%   hard_bits (N x 1 uint8)
%   success (logical)
%   iters (scalar)
%   Lpost_out (N x 1 double)
%   G (updated messages)

    if nargin < 4 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'max_iter'),   cfg.max_iter = 50; end
    if ~isfield(cfg,'saturation'), cfg.saturation = 25.0; end
    if ~isfield(cfg,'beta_code'),  cfg.beta_code = 1.0; end

    M = G.M; N = G.N;
    syndrome = uint8(syndrome(:));
    priors   = double(priors(:));

    if numel(syndrome) ~= M, error('Syndrome size mismatch.'); end
    if numel(priors)   ~= N, error('Priors size mismatch.'); end

    sat = cfg.saturation;

    % ----- initialize messages -----
    % Common choice: initialize v2c with priors, c2v with zeros.
    for v = 1:N
        val = saturate(priors(v), sat);
        edges_v = G.vn_edge{v};
        if isempty(edges_v), continue; end
        G.v2c(edges_v) = val;
        G.c2v(edges_v) = 0.0;
    end

    hard_bits = zeros(N,1,'uint8');
    Lpost_out = zeros(N,1);
    success = false;

    for iter = 1:cfg.max_iter

        % =========================
        % (1) CN update: all checks
        % =========================
        for c = 1:M
            edges = G.cn_edge{c};
            d = numel(edges);
            if d == 0, continue; end

            msgs = G.v2c(edges);
            abs_msgs = abs(msgs);

            if d == 1
                % Single-neighbor check: with no special rule, the outgoing
                % magnitude is min over empty set. A common convention is 0.
                % Sign still enforces syndrome.
                out_sign = (syndrome(c) == 0) * 2 - 1; % 0->+1, 1->-1
                G.c2v(edges(1)) = saturate(out_sign * 0.0, sat);
                continue;
            end

            % min1/min2
            [min1, min_idx] = min(abs_msgs);
            abs2 = abs_msgs;
            abs2(min_idx) = inf;
            min2 = min(abs2);

            % robust sign convention: negative -> -1, else +1 (zero -> +1)
            sgn = ones(size(msgs));
            sgn(msgs < 0) = -1;
            sign_prod = prod(sgn);

            % incorporate syndrome sign: (-1)^{syndrome(c)}
            if syndrome(c) ~= 0
                sign_prod = -sign_prod;
            end

            for k = 1:d
                if k == min_idx
                    mag = min2;
                else
                    mag = min1;
                end

                self_sign = 1;
                if msgs(k) < 0
                    self_sign = -1;
                end
                out_sign = sign_prod * self_sign;

                G.c2v(edges(k)) = saturate(out_sign * mag * cfg.beta_code, sat);
            end
        end

        % =========================
        % (2) VN update: all vars
        % =========================
        for v = 1:N
            edges_v = G.vn_edge{v};
            if isempty(edges_v), continue; end

            sum_c2v = sum(G.c2v(edges_v));
            Lpost = priors(v) + sum_c2v;

            for ee = 1:numel(edges_v)
                e = edges_v(ee);
                extrinsic = Lpost - G.c2v(e);
                G.v2c(e) = saturate(extrinsic, sat);
            end
        end

        % =========================
        % (3) Hard decision + check
        % =========================
        for v = 1:N
            edges_v = G.vn_edge{v};
            if isempty(edges_v)
                Lpost = priors(v);
            else
                Lpost = priors(v) + sum(G.c2v(edges_v));
            end
            Lpost_out(v) = Lpost;
            hard_bits(v) = uint8(Lpost < 0);
        end

        success = syndrome_check(G, hard_bits, syndrome);
        if success
            iters = iter;
            return;
        end
    end

    iters = cfg.max_iter;
end

function y = saturate(x, limit)
    y = min(max(x, -limit), limit);
end

function ok = syndrome_check(G, hard_bits, syndrome)
    M = G.M;
    ok = true;
    for c = 1:M
        vs = G.cn_to_vn{c};
        parity = uint8(0);
        for kk = 1:numel(vs)
            parity = bitxor(parity, hard_bits(vs(kk)));
        end
        if parity ~= syndrome(c)
            ok = false;
            return;
        end
    end
end
