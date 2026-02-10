clear;
load('simulation_data.mat');
load('non_isomorphic_TS.mat');
H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);

p = 1e-3;
L = log((1-p)/p);

wmax  = 5;   % max weight to enumerate (budget)

strength = [];
critical = [];
out_k = 0;
keys = repVN_by_class.keys;

for k = 1:length(repVN_by_class)
    supp = repVN_by_class(keys{k});
    subH = H(:, supp);
    rkeep = (sum(subH, 2) ~= 0);
    subH = subH(rkeep, :);

    subH = spones(subH) ~= 0;   % logical sparse

    dv = full(sum(subH, 1)).';       % n x 1
    dv_vec = full(subH * (dv - 1));  % M x 1

    llrs = L * ones(1, numel(supp));

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
        out_k = out_k + 1;
        TS{out_k} = supp;
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

nShow = length(TS);
for i = 1 : nShow
TSrep = TS{i};
plot_ts_from_H(H, TSrep, struct( ...
'show_vn_lbl',true, ...
'show_cn_lbl',true, ...
'unsat_def',"odd"));
end
