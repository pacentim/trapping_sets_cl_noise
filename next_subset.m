function [subset, last_generated_indicator] = next_subset(n, k, first_to_generate_indicator)
%NEXT_SUBSET Generate next k-subset of {1,...,n} (lex-like, via NEXKSB state).
% Returns:
%   subset : 1 x k vector of indices
%   last_generated_indicator : true iff this subset is the last one for this (n,k)

    persistent mem_old ind

    if nargin < 3, first_to_generate_indicator = 0; end

    % Reinitialize when called first time for this k-stream
    if isempty(mem_old) || first_to_generate_indicator == 1
        mem_old = initialize(n, k);
        ind = 0;
    end

    mem = NEXKSB(n, k, mem_old);

    mem_old = mem;
    ind = ind + 1;

    % "last generated" logic as in your code, but after state update
    last_generated_indicator = ~(mem_old.mtc == 1 || ind < 1);

    subset = mem.a(1:k);  % IMPORTANT: only first k entries are the subset

    if last_generated_indicator
        clear mem_old ind
    end
end


function mem = NEXKSB(n, k, mem_old)
%NEXKSB One step update of the k-subset state.
% State uses a(1:k) subset and a(k+1)=n+1 sentinel.

    a   = mem_old.a;
    mtc = mem_old.mtc;
    h   = mem_old.h;
    m2  = mem_old.m2;

    if mtc == 0
        m2o = 0;
        ho  = k;
    else
        if m2 < n - h
            ho = 0;
        else
            ho = h;
        end
        ho = ho + 1;

        % With sentinel a(k+1)=n+1, this index is always in 1..k
        m2o = a(k + 1 - ho);
    end

    ao = a;

    % Update the tail; indices are valid because ao has length k+1
    for j = 1:ho
        ao(k + j - ho) = m2o + j;
    end

    % Maintain sentinel explicitly
    ao(k+1) = n + 1;

    if ao(1) ~= n - k + 1
        mtco = 1;
    else
        mtco = 0;
    end

    mem.a   = ao;
    mem.mtc = mtco;
    mem.h   = ho;
    mem.m2  = m2o;
end


function mem = initialize(n, k)
%INITIALIZE Initialize the NEXKSB state with sentinel.
    a = [1:k, n+1];   % sentinel at position k+1
    mtc = 0;
    h = k;
    m2 = 0;

    mem.a   = a;
    mem.mtc = mtc;
    mem.h   = h;
    mem.m2  = m2;
end
