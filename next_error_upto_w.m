function [subset, wgt, done] = next_error_upto_w(n, wmax, reset)
%NEXT_ERROR_UPTO_W Enumerate supports of weight 0..wmax using next_subset.

    persistent cur_w stage finished

    if nargin < 3, reset = false; end

    if reset || isempty(cur_w)
        cur_w = 0;
        stage = 0;      % 0 -> return weight-0 once, then move to weight 1
        finished = false;
    end

    if finished
        subset = [];
        wgt = [];
        done = true;
        return;
    end

    % Weight 0 once
    if cur_w == 0 && stage == 0
        subset = [];
        wgt = 0;
        done = false;
        stage = 1;
        if wmax == 0
            finished = true;
        else
            cur_w = 1;
        end
        return;
    end

    % For cur_w >= 1: drive next_subset stream
    if stage == 1
        first = 1;  % initialize stream for this weight
        stage = 2;
    else
        first = 0;
    end

    [subset, last] = next_subset(n, cur_w, first);
    wgt = cur_w;
    done = false;

    if last
        cur_w = cur_w + 1;
        stage = 1; % next call initializes next weight stream
        if cur_w > wmax || cur_w==n
            finished = true;
        end
    end
end
