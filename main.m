clear;
load('simulation_data.mat');
H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);
% dv_avg = sum(sum(H,1))/990;

[NL, CL] = FindNeighbors(H);
Lp = enumeratePaths(NL,CL,5,5,5);
[~, cycles] = EnumerateCycles_Bani_NON_QC(H, 10, Lp);
cycles = uniquecell(cycles);

% Neighbor lists once
[NV, NC] = tanner_neighbors_from_H(H);

p = 1e-3;
L = log((1-p)/p);

kappa = 20;
wmax  = 5;   % max weight to enumerate (budget)

TS = cell(1, numel(cycles));
strength = nan(1, numel(cycles));
critical = nan(1, numel(cycles));
out_k = 0;

for k = 1:numel(cycles)
    supp = cycles{k};
    exp_supp = expand_one_cycle(NV, NC, supp, kappa);

    subH = H(:, exp_supp);
    rkeep = (sum(subH, 2) ~= 0);
    subH = subH(rkeep, :);

    subH = spones(subH) ~= 0;   % logical sparse

    dv = full(sum(subH, 1)).';       % n x 1
    dv_vec = full(subH * (dv - 1));  % M x 1

    llrs = L * ones(1, numel(exp_supp));

    M = size(subH,1);
    n = size(subH,2);

    % Union bitset of failing indices (local 1..n)
    nWords = ceil(n/64);
    union_idx = false(1,n);

    minw = inf;
    ncrit = 0;
    fail_count_by_w = zeros(1, wmax+1);  % weights 0..wmax

    % Syndrome for current subset (recomputed each time; safe + simple)
    syn = false(M,1);

    % Enumerate subsets up to weight wmax
    reset = true;
    while true
        [S, w, done] = next_error_upto_w(n, wmax, reset);
        reset = false;
        if done, break; end

        % Build syndrome for this support S
        syn(:) = false;
        for t = 1:numel(S)
            syn = xor(syn, subH(:, S(t)));
        end

        syndrome = double(syn).';
        [~, success] = expcontr_minsum_decode_sparse(subH, syndrome, llrs, [], dv_vec);

        if ~success
            union_idx(S) = (union_idx(S)+1)>0;
            fail_count_by_w(w + 1) = fail_count_by_w(w + 1) + 1;

            if w < minw
                minw = w;
                ncrit = 1;
            elseif w == minw
                ncrit = ncrit + 1;
            end
        end
    end

    if isfinite(minw)
        TS_this = exp_supp(union_idx);

        out_k = out_k + 1;
        TS{out_k} = TS_this;
        strength(out_k) = minw;
        critical(out_k) = ncrit;

        % Optional:
        % fprintf("k=%d: n=%d, wmax=%d, minw=%d, ncrit=%d, fail_w1=%d\n", ...
        %     k, n, wmax, minw, ncrit, fail_count_by_w(1+1));
    end
end

TS = TS(1:out_k);
strength = strength(1:out_k);
critical = critical(1:out_k);

% Canonicalize each TS (sorted row vectors)
TS_sorted = cellfun(@(x) sort(x(:).'), TS, 'UniformOutput', false);

% Build unique keys (string representation)
keys = cellfun(@(x) sprintf('%d,', x), TS_sorted, 'UniformOutput', false);

% Unique with stable order
[~, ia] = unique(keys, 'stable');

% Keep only unique entries
TS = TS(ia);
strength = strength(ia);
critical = critical(ia);
