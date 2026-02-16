function rows = generate_monotone_rows(levels)
% All length-7 nondecreasing sequences over given sorted levels.
% Here levels = [-3..3], so count is C(13,7)=1716.

    L = numel(levels);  % 7
    n = 7;

    rows = zeros(0,n,'int8');
    idx = ones(1,n); % indices into levels

    while true
        rows(end+1,:) = levels(idx); %#ok<AGROW>

        % next nondecreasing multiset combination
        p = n;
        while p >= 1 && idx(p) == L
            p = p - 1;
        end
        if p == 0
            break;
        end
        idx(p) = idx(p) + 1;
        idx(p+1:end) = idx(p);
    end
end