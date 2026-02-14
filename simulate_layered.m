%% Test layered decoder

clear;
load('simulation_data.mat');
load('failure_sets.mat');
H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);
T = build_tanner_graph(H);
wmax = 5;
p = 1e-3;
L = log((1-p)/p);
layers = split_graph_into_layers(T,100);
corrects = [];
llrs = [L * ones(1,990) zeros(1,size(H,2)-990)];

for i = 1 : length(TS)
    supp = TS{i};
    n = length(supp);
    subH = H(:,supp);
     % Union bitset of failing indices (local 1..n)
    nWords = ceil(n/64);
    union_idx = false(1,n);

    minw = inf;
    ncrit = 0;
    fail_count_by_w = zeros(1, wmax+1);  % weights 0..wmax

    % Syndrome for current subset (recomputed each time; safe + simple)
    syn = false(T.M,1);

    % Enumerate subsets up to weight wmax
    reset = true;
    while true
        [S, w, done] = next_error_upto_w(n, wmax, reset);
        reset = false;
        if done, break; end

        % Build syndrome for this support S
        syn(:) = false;
        for t = 1:numel(S)
            % syn = xor(syn, fullsubH(:, S(t)));
            syn = xor(syn, subH(:, S(t)));
        end

        syndrome = double(syn).';
        % [~, success] = super_layered_decode(T, layers, syndrome, llrs);
        [~, success] = super_layered_decode_faid(T, layers, syndrome, llrs);
        if ~success
            break;
        end

    end

    if success
        corrects = [corrects 1];
    else
        corrects = [corrects 0];
    end
end

