function [N_v, N_c] = tanner_neighbors_from_H(H, do_check)
% [N_v, N_c] = tanner_neighbors_from_H(H, do_check)
% Build Tanner-graph neighbor lists from binary parity-check matrix H.
%
% H is m-by-n with entries in {0,1} (logical or numeric).
% - N_v{v} = list of check-node indices (row indices) connected to VN v (column v)
% - N_c{c} = list of variable-node indices (col indices) connected to CN c (row c)
%
% do_check (optional, default true): verify N_v and N_c are consistent.

    if nargin < 2
        do_check = true;
    end

    if ~ismatrix(H)
        error('H must be a 2-D matrix.');
    end

    % Ensure logical for fast find
    H = (H ~= 0);

    [m, n] = size(H);

    % Preallocate cell arrays
    N_v = cell(1, n);
    N_c = cell(1, m);

    % CN neighbors: for each check (row), list variable indices
    for c = 1:m
        N_c{c} = find(H(c, :));   % row vector of VN indices
    end

    % VN neighbors: for each variable (col), list check indices
    for v = 1:n
        N_v{v} = find(H(:, v)).';  % row vector of CN indices
    end

    if do_check
        % Consistency: c in N_v{v} iff v in N_c{c}
        for v = 1:n
            for c = N_v{v}
                if ~any(N_c{c} == v)
                    error('Inconsistency: CN %d in N_v{%d} but VN %d not in N_c{%d}.', c, v, v, c);
                end
            end
        end
        for c = 1:m
            for v = N_c{c}
                if ~any(N_v{v} == c)
                    error('Inconsistency: VN %d in N_c{%d} but CN %d not in N_v{%d}.', v, c, c, v);
                end
            end
        end
    end
end
