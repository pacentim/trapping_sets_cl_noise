function [Vcore, stats] = peel_dangling_vns(N_v, N_c, Vinit)
%PEEL_DANGLING_VNS Iteratively peel VNs incident to any induced degree-1 CN.
%
% Inputs:
%   N_v   : cell array, N_v{v} = CN neighbors of VN v (global indexing)
%   N_c   : cell array, N_c{c} = VN neighbors of CN c (global indexing)
%   Vinit : row/col vector of VN indices (global)
%
% Outputs:
%   Vcore : remaining VN set after peeling (possibly empty)
%   stats : struct with fields rounds, removed_per_round

    nVar = numel(N_v);
    nChk = numel(N_c);

    V = unique(Vinit(:)).';
    if isempty(V)
        Vcore = V;
        stats.rounds = 0;
        stats.removed_per_round = [];
        return;
    end

    Vmark = zeros(1, nVar, 'uint32');
    Cmark = zeros(1, nChk, 'uint32');
    tV = uint32(1);
    tC = uint32(1);

    removed = zeros(1, 64, 'uint32');  % will grow if needed
    r = 0;

    while true
        % Mark current V
        tV = tV + 1;
        Vmark(V) = tV;

        % CN neighborhood of current V
        C = unique([N_v{V}]);
        if isempty(C)
            break; % no checks, nothing dangling
        end

        % Induced degrees of checks w.r.t. V
        deg_ind = zeros(1, numel(C), 'uint16');
        for ii = 1:numel(C)
            c = C(ii);
            deg_ind(ii) = sum(Vmark(N_c{c}) == tV);
        end
        C1 = C(deg_ind == 1); % dangling checks

        if isempty(C1)
            break; % no dangling checks -> fixed point
        end

        % Mark dangling checks
        tC = tC + 1;
        Cmark(C1) = tC;

        % Peel: remove VNs having >=1 neighbor in C1
        peelMask = false(1, numel(V));
        for i = 1:numel(V)
            v = V(i);
            peelMask(i) = any(Cmark(N_v{v}) == tC);
        end

        if ~any(peelMask)
            break; % should not happen, but safe
        end

        r = r + 1;
        if r > numel(removed)
            removed = [removed, zeros(1, numel(removed), 'uint32')]; %#ok<AGROW>
        end
        removed(r) = uint32(nnz(peelMask));

        V = V(~peelMask);
        if isempty(V)
            break; % TS emptied by peeling
        end
    end

    Vcore = V;
    stats.rounds = r;
    stats.removed_per_round = double(removed(1:r));
end
