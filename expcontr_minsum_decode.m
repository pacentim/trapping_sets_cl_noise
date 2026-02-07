function [hard_bits, success, iters, Lpost_out, G] = expcontr_minsum_decode(G, syndrome, priors, cfg, dv_ext)
% Flooding min-sum with degree-1 CN injection via recursion theta^(l).
%
% dv_ext:
%   - scalar: use same external VN degree for every degree-1 CN
%   - (M x 1) vector: dv_ext(c) used for check c

    if nargin < 4 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'max_iter'),   cfg.max_iter = 50; end
    if ~isfield(cfg,'saturation'), cfg.saturation = 25.0; end
    if ~isfield(cfg,'beta_code'),  cfg.beta_code = 1.0; end

    if nargin < 5 || isempty(dv_ext)
        error('Provide dv_ext (scalar or length-M vector) for degree-1 CN injection.');
    end

    M = G.M; N = G.N;
    syndrome = syndrome(:);
    priors   = priors(:);

    if numel(syndrome) ~= M, error('Syndrome size mismatch.'); end
    if numel(priors) ~= N,   error('Priors size mismatch.'); end

    sat = cfg.saturation;
    Y = priors(1);

    % normalize dv_ext
    if isscalar(dv_ext)
        dv_mode = 'scalar';
        dv_scalar = double(dv_ext);
    else
        dv_mode = 'vector';
        dv_vec = double(dv_ext(:));
        if numel(dv_vec) ~= M
            error('dv_ext vector must have length M.');
        end
    end

    % ----- initialize messages -----
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

        % =========================
        % (0) Update injected theta^(iter-1) -> theta^(iter)
        % =========================
        % Define theta^(0)=Y. After that: theta^(l)=sat(Y+(dv-1)theta^(l-1)).
        if iter == 1
            if strcmp(dv_mode,'scalar')
                theta_iter_scalar = saturate(Y, sat); % theta^(0)
            else
                theta_iter = saturate(Y * ones(M,1), sat);
            end
        else
            if strcmp(dv_mode,'scalar')
                theta_iter_scalar = saturate(Y + (dv_scalar-1) * theta_iter_scalar, sat);
            else
                theta_iter = saturate(Y + (dv_vec-1) .* theta_iter, sat);
            end
        end

        % =========================
        % (1) CN update: all checks
        % =========================
        for c = 1:M
            edges = G.cn_edge{c};
            d = numel(edges);
            if d == 0, continue; end

            if d == 1
                % ---- degree-1 CN: outgoing uses ONLY the injected external message ----
                if strcmp(dv_mode,'scalar')
                    mag = theta_iter_scalar;
                else
                    mag = theta_iter(c);
                end

                out_sign = 1;
                if syndrome(c) ~= 0
                    out_sign = -1;
                end

                G.c2v(edges(1)) = saturate(out_sign * mag * cfg.beta_code, sat);
                continue;
            end

            % ---- standard min-sum for d >= 2 ----
            msgs = G.v2c(edges);
            abs_msgs = abs(msgs);

            [min1, min_idx] = min(abs_msgs);
            abs2 = abs_msgs;
            abs2(min_idx) = inf;
            min2 = min(abs2);

            sgn = ones(size(msgs));
            sgn(msgs < 0) = -1;
            sign_prod = prod(sgn);

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
                if G.v2c(edges(k)) < 0
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
            Lpost = priors(v) + sum(G.c2v(edges_v));
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
        if parity ~= uint8(syndrome(c))
            ok = false;
            return;
        end
    end
end
