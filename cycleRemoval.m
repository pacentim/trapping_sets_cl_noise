%% Removing 4-cycles: 
% S. Sankaranarayanan and B. Vasic, "Iterative decoding of linear block codes: a parity-check orthogonalization approach," in IEEE Transactions on Information Theory, vol. 51, no. 9, pp. 3347-3353, Sept. 2005
function Hnew = cycleRemoval(H)
% REMOVE4CYCLES  Eliminate all 4-cycles using the auxiliary-variable/check
%                transformation described in your screenshot.
%
%   Inputs:
%     H : m-by-n binary parity-check matrix (checks x variables). Prefer sparse.
%
%   Outputs:
%     Hnew        : transformed matrix with no 4-cycles
%     info        : struct with iteration log
%
%   Notes:
%     - Implements Steps 1–4:
%         * pick a variable pair (u,v) with >=2 shared checks,
%         * wv := columnwise AND of columns u and v,
%         * add B that XORs wv into columns u and v (removes their shared
%           1s at the common checks),
%         * append a new variable column equal to wv,
%         * append a new check row S' with 1s at positions u, v, and the new column.
%
%     - All arithmetic is over GF(2).

    % Ensure binary & sparse
    if ~issparse(H), H = sparse(H ~= 0); else, H = spones(H); end

    l = 0;
    H_l = H;

    itPairs   = [];   % [u v] chosen at each iteration
    itSharedK = [];   % how many shared checks before fixing (k = (H'*H)(u,v))
    itBroken  = [];   % how many 4-cycles broken in that step

    while true
        % Step 1: compute H_l' * H_l and find pairs with >=2 shared checks
        S = H_l' * H_l;                % variable-variable length-2 paths
        S = S - diag(diag(S));         % zero diagonal
        [uIdx, vIdx, kVals] = find(triu(S, 1));

        if isempty(kVals) || max(kVals) < 2
            % four-cycle free -> Step 4
            Hnew = H_l;
            break
        end

        % --- Choose a pair (u,v) involved in at least one 4-cycle ---
        % Heuristic: pick the pair with largest k to break many cycles at once.
        [kmax, pos] = max(kVals);
        u = uIdx(pos);  v = vIdx(pos);

        % Step 2: construct wv = u .* v (componentwise product over GF(2))
        wv = H_l(:,u) & H_l(:,v);      % logical AND == componentwise product

        % Number of cycles broken for this pair equals nchoosek(kmax,2)
        brokenHere = kmax*(kmax-1)/2;

        % Build B that has wv in columns u and v, zeros elsewhere
        [m_l, n_l] = size(H_l);
        B = sparse(m_l, n_l);
        B(:,u) = wv;
        B(:,v) = wv;

        % Step 3: apply transform
        % H_{l+1} upper-left block: H_l + B over GF(2)
        H_next = xor_sparse(H_l, B);

        % Append the new variable column (n_l+1) equal to wv
        H_next = [H_next, wv];

        % Append the auxiliary check S' with 1s at u, v, and (n_l+1)
        Sprime = sparse(1, n_l+1);
        Sprime(1, u) = 1;
        Sprime(1, v) = 1;
        Sprime(1, n_l+1) = 1;
        H_next = [H_next; Sprime];

        % Bookkeeping
        l = l + 1;
        itPairs(end+1,:)   = [u v];         %#ok<AGROW>
        itSharedK(end+1,1) = full(kmax);    %#ok<AGROW>
        itBroken(end+1,1)  = brokenHere;    %#ok<AGROW>

        % Continue
        H_l = spones(H_next);
    end

    % Report
    % info.iters   = l;
    % info.pairs   = itPairs;
    % info.kvals   = itSharedK;
    % info.broken  = itBroken;
    % info.final4c = count4cycles(Hnew);
    Hnew=full(Hnew);
end

% --------- helpers ---------

function C = xor_sparse(A, B)
% XOR for sparse 0/1 matrices (GF(2) addition).
    C = A + B;
    % In GF(2), values 2 -> 0, values 1 -> 1.
    % Efficient way: mod 2 and sparsify.
    [i,j,s] = find(C);
    s = mod(s,2);
    keep = s ~= 0;
    C = sparse(i(keep), j(keep), s(keep), size(A,1), size(A,2));
end
