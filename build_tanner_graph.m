function G = build_tanner_graph(H)
%BUILD_TANNER_GRAPH Build Tanner graph structures from a binary parity-check matrix H.
% Mimics the C++ TannerGraph build_from() layout and edge indexing.
%
% Input:
%   H : (M x N) logical/double/sparse with entries in {0,1}
%
% Output struct G fields:
%   M, N, E
%   cn_to_vn : cell(M,1), each is row vector of VN indices (1-based)
%   vn_to_cn : cell(N,1), each is row vector of CN indices (1-based)
%   cn_edge  : cell(M,1), each is row vector of edge ids (1..E), aligned to cn_to_vn
%   vn_edge  : cell(N,1), each is row vector of edge ids (1..E), aligned to vn_to_cn
%   v2c, c2v  : (E x 1) double messages

    if ~issparse(H)
        % sparse is typically faster for neighborhood enumeration
        H = sparse(H);
    end

    [M,N] = size(H);
    G.M = M;
    G.N = N;

    % Neighbors
    G.cn_to_vn = cell(M,1);
    G.vn_to_cn = cell(N,1);

    % First pass: collect neighbors and count edges
    edge_count = 0;
    for i = 1:M
        js = find(H(i,:));            % VN indices in row i
        G.cn_to_vn{i} = js(:).';      % row vector
        edge_count = edge_count + numel(js);
        for kk = 1:numel(js)
            j = js(kk);
            G.vn_to_cn{j} = [G.vn_to_cn{j}, i];
        end
    end
    G.E = edge_count;

    % Allocate edge index maps aligned with neighbor lists
    G.cn_edge = cell(M,1);
    for i = 1:M
        G.cn_edge{i} = zeros(1, numel(G.cn_to_vn{i}), 'int32');
    end

    G.vn_edge = cell(N,1);
    for j = 1:N
        G.vn_edge{j} = zeros(1, numel(G.vn_to_cn{j}), 'int32');
    end

    % Second pass: assign global edge ids in the same order as C++:
    % loop CN i=1..M, then its neighbor list order, edge++.
    e = 1;
    for i = 1:M
        js = G.cn_to_vn{i};
        for k = 1:numel(js)
            j = js(k);

            G.cn_edge{i}(k) = int32(e);

            % find position of CN i inside vn_to_cn{j}
            cn_list = G.vn_to_cn{j};
            pos = find(cn_list == i, 1, 'first');
            if isempty(pos)
                error('TannerGraph build_from: graph inconsistency at CN %d, VN %d', i, j);
            end
            G.vn_edge{j}(pos) = int32(e);

            e = e + 1;
        end
    end

    if e-1 ~= G.E
        error('TannerGraph build_from: edge count mismatch');
    end

    % Messages
    G.v2c = zeros(G.E,1);
    G.c2v = zeros(G.E,1);
end
