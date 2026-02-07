function Vexp = expand_one_cycle(N_v, N_c, V0, kappa)
%EXPAND_ONE_CYCLE_FAST_TS Timestamp-marking version: avoids clearing masks.

    nVar = numel(N_v);
    nChk = numel(N_c);

    V = unique(V0(:)).';

    Vmark = zeros(1, nVar, 'uint32');
    Cmark = zeros(1, nChk, 'uint32');
    tV = uint32(1);
    tC = uint32(1);

    while true
        % Mark V
        Vmark(V) = tV;

        % C = neighbors of V
        C = unique([N_v{V}]);
        Cmark(C) = tC;

        % Induced degrees
        deg_ind = zeros(1, numel(C), 'uint16');
        for ii = 1:numel(C)
            c = C(ii);
            deg_ind(ii) = sum(Vmark(N_c{c}) == tV);
        end
        C1 = C(deg_ind == 1);

        % Collect W
        Wcells = cell(1, numel(C1));
        widx = 0;

        for c = C1
            neighV = N_c{c};
            cand = neighV(Vmark(neighV) ~= tV);  % not in V
            if isempty(cand), continue; end

            keep = false(1, numel(cand));
            for jj = 1:numel(cand)
                w = cand(jj);
                cn_w = N_v{w};
                dv = length(cn_w);

                % condition: at least two checks in C besides c
                % <=> total checks in C among cn_w is at least 3 (includes c)
                if sum(Cmark(cn_w) == tC) >= 2 %(floor(dv/2)+1)
                    keep(jj) = true;
                end
            end

            if any(keep)
                widx = widx + 1;
                Wcells{widx} = cand(keep);
            end
        end

        if widx == 0
            break;
        end

        W = unique([Wcells{1:widx}]);
        if isempty(W)
            break;
        end

        Vnew = unique([V, W]);
        if numel(Vnew) == numel(V)
            break;
        end
        V = Vnew;

        if numel(V) >= kappa
            
            break;
        end

        % advance timestamps (no clearing)
        tV = tV + 1;
        tC = tC + 1;
    end

    Vexp = V;
end
