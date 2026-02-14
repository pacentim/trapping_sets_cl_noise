function layers = split_graph_into_layers(G, max_per_layer)
%SPLIT_GRAPH_INTO_LAYERS MATLAB port of the provided Python function.
%
% layers = split_graph_into_layers(G, max_per_layer)
%
% Inputs:
%   G: Tanner graph struct with fields:
%      G.M, G.cn_to_vn (cell array, cn_to_vn{c} is vector of VN indices, 1-based)
%   max_per_layer: positive integer
%
% Output:
%   layers: cell array, layers{ell} = vector of CN indices (1-based)

    if nargin < 2 || isempty(max_per_layer) || max_per_layer < 1 || floor(max_per_layer) ~= max_per_layer
        error("max_per_layer must be a positive integer.");
    end

    M = double(G.M);

    % ---- cn_vars and degrees ----
    cn_vars = cell(M,1);
    deg = zeros(M,1);

    for c = 1:M
        neigh = G.cn_to_vn{c};
        if isempty(neigh)
            cn_vars{c} = zeros(1,0);
            deg(c) = 0;
        else
            arr = unique(double(neigh(:)).');   % unique VNs (row vector)
            cn_vars{c} = arr;
            deg(c) = numel(arr);
        end
    end

    % ---- sort by (-deg, +index) like np.lexsort((order, -deg)) ----
    order = (1:M).';
    T = table(-deg, order);              % primary: -deg ascending (=> deg descending), secondary: order ascending
    T = sortrows(T, {'Var1','order'});   % Var1 is -deg
    check_order = T.order;               % column vector of CN ids (1..M)

    % ---- greedy layering with max_per_layer and best-fit (min size) ----
    layers = {};             % cell array of vectors
    layer_used = {};         % cell array of logical masks (N x 1) for fast membership
    layer_sizes = [];        % numeric vector

    N = double(G.N);

    for idx = 1:numel(check_order)
        c = check_order(idx);
        vars_c = cn_vars{c};

        best_ell = 0;        % 0 means "not found" (MATLAB)
        best_size = inf;

        for ell = 1:numel(layers)
            if layer_sizes(ell) >= max_per_layer
                continue;
            end

            if isempty(vars_c)
                feasible = true;
            else
                used_mask = layer_used{ell};
                feasible = ~any(used_mask(vars_c));
            end

            if feasible
                if layer_sizes(ell) < best_size
                    best_size = layer_sizes(ell);
                    best_ell = ell;
                end
            end
        end

        if best_ell == 0
            % create new layer
            layers{end+1} = c; %#ok<AGROW>
            used_mask = false(N,1);
            if ~isempty(vars_c)
                used_mask(vars_c) = true;
            end
            layer_used{end+1} = used_mask; %#ok<AGROW>
            layer_sizes(end+1) = 1; %#ok<AGROW>
        else
            layers{best_ell}(end+1) = c;
            if ~isempty(vars_c)
                layer_used{best_ell}(vars_c) = true;
            end
            layer_sizes(best_ell) = layer_sizes(best_ell) + 1;
        end
    end

    % store as row vectors (cosmetic)
    for ell = 1:numel(layers)
        layers{ell} = layers{ell}(:).';
    end
end
