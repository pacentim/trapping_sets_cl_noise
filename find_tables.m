%% Exhaustive LUT2 search for ONE fixed TS (dv=2), keep LUTs that correct ALL errors up to wmax

clear;
load('simulation_data.mat');
load('failure_sets.mat');

H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);

T = build_tanner_graph(H);
layers = split_graph_into_layers(T, 100);

wmax = 3;
p = 1e-3;
Lch = log((1-p)/p);
llrs = [Lch * ones(1,990) zeros(1,size(H,2)-990)];

% Enumerate all monotone LUT2neg rows (1716 candidates)
levels = int8(-3:3);
rows = generate_monotone_rows(levels);   % 1716 x 7
% rows = [[-3 -3 -3 -3 -3 -1 2];
%  [-3 -3 -3 -3 -3 1 2];
% [-3 -3 -3 -3 -2 1 2];
% [-3 -3 -3 -3 0 1 2];
% [-3 -3 -3 -2 -2 1 2];
% [-3 -3 -3 -2 -1 1 3];
% [-3 -3 -2 -2 0 2 2]];

% -----------------------
% Choose ONE trapping set
% -----------------------
ts_idx = 4;              % <<< set this
supp = TS{ts_idx};
n = numel(supp);
subH = H(:, supp);       % M x n (logical/double)

% Pre-allocate result container
good_rows = zeros(0,7,'int8');
good_k    = [];          % indices in "rows"

% For syndrome computation
syn = false(T.M,1);

% Optional: default LUT3 (leave as your default inside decoder)
% We'll pass only LUT2neg via 'lut' struct
lut = struct();
cfg = struct(); % optionally set cfg.max_iter, cfg.prior_scale, cfg.layer_damping, etc.

fprintf('TS #%d: n=%d, testing %d monotone LUT2 rows, wmax=%d\n', ts_idx, n, size(rows,1), wmax);

for r = 1:size(rows,1)
    LUT2neg = rows(r,:);

    % set the LUT for dv=2
    lut.LUT2neg = LUT2neg;

    % Assume good until proven otherwise
    is_good = true;

    % Enumerate all error supports S up to weight wmax
    reset = true;
    while true
        [S, w, done] = next_error_upto_w(n, wmax, reset);
        reset = false;
        if done
            break;
        end

        % Build syndrome for this support S (in FULL graph coordinates)
        syn(:) = false;
        for t = 1:numel(S)
            syn = xor(syn, subH(:, S(t)));
        end

        syndrome = double(syn).';

        % Decode using FAID with this LUT
        % IMPORTANT: super_layered_decode_faid must accept (cfg,state,lut) or (cfg,state,lutStruct)
        % If your modified signature is (..., cfg, state, lut), use:
        [~, success] = super_layered_decode_faid(T, layers, syndrome, llrs, cfg, [], lut);

        if ~success
            is_good = false;
            break; % discard this LUT immediately
        end
    end

    if is_good
        good_rows(end+1,:) = LUT2neg; %#ok<AGROW>
        good_k(end+1) = r;            %#ok<AGROW>
        fprintf('GOOD LUT found (r=%d): %s\n', r, mat2str(double(LUT2neg)));
    end

    % Progress
    if mod(r, 100) == 0
        fprintf('Tested %d / %d, good so far: %d\n', r, size(rows,1), size(good_rows,1));
    end
end

fprintf('Done. Good LUTs: %d / %d\n', size(good_rows,1), size(rows,1));

% Save results
results = struct();
results.ts_idx = ts_idx;
results.supp = supp;
results.wmax = wmax;
results.p = p;
results.Lch = Lch;
results.good_rows = good_rows;
results.good_k = good_k;
save('good_lut2_rows_for_ts.mat', 'results');

%% ---------------- helper: monotone row generator ----------------
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
