function sub = build_subgraph(N_v, N_c, Vset)
    V = sort(unique(Vset(:)))';
    C = unique([N_v{V}]); C = sort(C(:))';

    % maps: global id -> local id
    Vmax = max(V); Cmax = max(C);
    Vmap = zeros(1, Vmax); Vmap(V) = 1:numel(V);
    Cmap = zeros(1, Cmax); Cmap(C) = 1:numel(C);

    % local adjacency lists
    Nv_local = cell(1, numel(V));
    for j = 1:numel(V)
        v = V(j);
        cn = N_v{v};
        cn = cn(ismember(cn, C));
        Nv_local{j} = Cmap(cn);
    end

    Nc_local = cell(1, numel(C));
    for i = 1:numel(C)
        c = C(i);
        vn = N_c{c};
        vn = vn(ismember(vn, V));
        Nc_local{i} = Vmap(vn);
    end

    sub.V = V;
    sub.C = C;
    sub.Nv = Nv_local;
    sub.Nc = Nc_local;
end

function sub = annotate_degree1_checks(sub)
    deg_ind = zeros(1, numel(sub.C));
    for i = 1:numel(sub.C)
        deg_ind(i) = numel(sub.Nc{i}); % induced degree inside subgraph
    end
    sub.degC_ind = deg_ind;
    sub.C1_local = find(deg_ind == 1); % local indices of degree-1 CNs
end
