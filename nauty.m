
clear;
load('simulation_data.mat');
H = H(:,1:990);
H(sum(H,2)==0,:) = [];
H = cycleRemoval(H);
load('trapping_sets.mat');
outDir = "C:\temp";
if ~isfolder(outDir)
    mkdir(outDir);
end
out_g6 = "C:\temp\trapping_sets.g6";
% I = I(1:15);  % your trapping sets (cell array of VN index vectors)

[repVN_by_class, count_by_class, members_by_class, meta] = ...
    classify_ts_store_repVN(H, I, "/home/michele/nauty2_9_3", false);

fprintf("Total non-isomorphic classes: %d\n", repVN_by_class.Count);

keys = repVN_by_class.keys;
for i = 1:min(10, numel(keys))
    key = keys{i};
    rep = repVN_by_class(key);
    cnt = count_by_class(key);
    fprintf("Class %d: multiplicity=%d, rep VN = [%s]\n", ...
        i, cnt, num2str(rep));
end
% 
keys = repVN_by_class.keys;

nShow = min(20, numel(keys));
for i = 1 : nShow
    key   = keys{i};
    TSrep = repVN_by_class(key);

    fprintf("Class %d/%d: multiplicity=%d, |TS|=%d\n", ...
        i, numel(keys), count_by_class(key), numel(TSrep));

    plot_ts_from_H(H, TSrep, struct( ...
        'show_vn_lbl',true, ...
        'show_cn_lbl',true, ...
        'unsat_def',"odd"));
end
% % 
% % TSa = repVN_by_class(keys{11});
% % TSb = repVN_by_class(keys{13});
% % 
% [isIso, kA, kB, info] = are_iso_ts_labelg(H, TSa, TSb, "/home/michele/nauty2_9_3");
% 
% disp(info)
% fprintf("isIso = %d\n", isIso);
% fprintf("keyA = %s\nkeyB = %s\n", kA, kB);
% 
% inv1 = ts_invariants(H, TSa);
% inv2 = ts_invariants(H, TSb);
% 
% isequal(inv1.dv, inv2.dv)
% isequal(inv1.dc, inv2.dc)
% 
% 
% A1 = graph6_to_adj(kA);
% A2 = graph6_to_adj(kB);
% e1 = sort(eig(double(A1)));
% e2 = sort(eig(double(A2)));
% fprintf("Spectra equal? %d\n", max(abs(e1-e2)) < 1e-9);



% A = graph62adj(g6list(1));
% k = 5;                       % noto dal batch
% [A_vncn, k, m] = split_tanner_adj(A, k);






function p = win2wsl_path(pwin)
% Convert Windows path like C:\temp\file.g6 to /mnt/c/temp/file.g6
    pwin = char(pwin);
    drive = lower(pwin(1));
    rest  = strrep(pwin(3:end), '\', '/');
    p = sprintf('/mnt/%c/%s', drive, rest);
end





function [repVN_by_class, count_by_class, members_by_class, meta] = ...
    classify_ts_store_repVN(H, TSlist, nautyDirLinux, store_members)
%CLASSIFY_TS_STORE_REPVN
% Canonicalize induced Tanner subgraphs of trapping sets using nauty labelg
% with a VN/CN color constraint, and store one representative VN-index set
% per non-isomorphism class.
%
% Key points (robust):
%   - Group by (k,m), where k = |TS| and m = #CN touching TS (in induced subgraph)
%   - NO padding with isolated CN nodes
%   - Strip output lines from labelg; silence stderr to avoid key pollution
%
% Inputs:
%   H             : parity-check matrix (logical or 0/1)
%   TSlist        : cell array, TSlist{i} = vector of VN indices (1-based)
%   nautyDirLinux : e.g. "/home/michele/nauty2_9_3"
%   store_members : true/false; if true, store all VN sets per class
%
% Outputs:
%   repVN_by_class   : containers.Map key=canonical g6, value=VN index row vector
%   count_by_class   : containers.Map key=canonical g6, value=double count
%   members_by_class : containers.Map key=canonical g6, value=cell array of VN sets
%                      (empty [] if store_members=false)
%   meta             : struct with grouping + file info

    if nargin < 3 || strlength(nautyDirLinux) == 0
        nautyDirLinux = "/home/michele/nauty2_9_3";
    end
    if nargin < 4
        store_members = false;
    end

    labelg = nautyDirLinux + "/labelg";

    % Normalize H
    if ~islogical(H), H = (H ~= 0); end

    % Basic sanity on TSlist
    if ~iscell(TSlist), error("TSlist must be a cell array"); end
    nVar = size(H,2);
    for i = 1:numel(TSlist)
        if isempty(TSlist{i}), continue; end
        if any(TSlist{i} < 1) || any(TSlist{i} > nVar)
            error("TSlist{%d} has VN indices out of range [1..%d]. Check H vs TS indexing.", i, nVar);
        end
    end

    % Precompute (k,m) for each TS
    N = numel(TSlist);
    ks = zeros(N,1);
    ms = zeros(N,1);

    for i = 1:N
        TS = unique(TSlist{i}(:).');
        ks(i) = numel(TS);

        subH = H(:, TS);
        subH(sum(subH,2)==0,:) = [];
        ms(i) = size(subH,1);
    end

    pairs = unique([ks ms], 'rows');   % each row: [k m]

    repVN_by_class = containers.Map('KeyType','char','ValueType','any');
    count_by_class = containers.Map('KeyType','char','ValueType','double');

    if store_members
        members_by_class = containers.Map('KeyType','char','ValueType','any');
    else
        members_by_class = [];
    end

    meta = struct();
    meta.pairs = pairs;
    meta.files = struct();
    meta.nautyDirLinux = nautyDirLinux;

    outDir = "C:\temp\nauty_tmp";
    if ~isfolder(outDir), mkdir(outDir); end

    if ~isfolder(outDir), mkdir(outDir); end
    meta.outDir = outDir;

    for ii = 1:size(pairs,1)
        k = pairs(ii,1);
        m = pairs(ii,2);

        idx = find(ks==k & ms==m);

        if isempty(idx)
            continue;
        end

        % Write all graphs for this (k,m) into a .g6 file
        in_g6 = fullfile(outDir, sprintf("ts_k%d_m%d.g6", k, m));
        fid = fopen(in_g6, "w");
        if fid < 0, error("Cannot open %s", in_g6); end

        for t = 1:numel(idx)
            TS = unique(TSlist{idx(t)}(:).');
            A = ts2adj_tanner_nopad(H, TS);    % n = k+m
            g6 = adj2graph6(A);
            fprintf(fid, "%s\n", g6);
        end
        fclose(fid);

        % Canonicalize with VN/CN partition constraint: first k are 'a', last m are 'b'
        fopt = sprintf('-fa^%db^%d', k, m);

        % IMPORTANT:
        % - Quote paths
        % - Strip stderr to avoid junk lines mixing into stdout
        in_g6_wsl = win2wsl_path(in_g6);

        cmd = sprintf('wsl bash -lc "%s -q %s \\"%s\\" 2>/dev/null"', ...
        labelg, fopt, in_g6_wsl);


        [status,out] = system(cmd);
        if status ~= 0
            error("labelg failed for k=%d,m=%d:\n%s", k, m, out);
        end

        canon = string(splitlines(string(out)));
        canon = strip(canon);
        canon = canon(canon ~= "");

        if numel(canon) ~= numel(idx)
            error("labelg output lines (%d) != input TS count (%d) for k=%d,m=%d", ...
                numel(canon), numel(idx), k, m);
        end

        % Update maps
        for t = 1:numel(idx)
            key = char(canon(t));              % canonical g6
            TS  = unique(TSlist{idx(t)}(:).'); % VN indices (original H indexing)

            if ~isKey(repVN_by_class, key)
                repVN_by_class(key) = TS;
                count_by_class(key) = 1;

                if store_members
                    members_by_class(key) = {TS};
                end
            else
                count_by_class(key) = count_by_class(key) + 1;

                if store_members
                    L = members_by_class(key);
                    L{end+1} = TS;
                    members_by_class(key) = L;
                end
            end
        end

        tag = sprintf("k%d_m%d", k, m);
        meta.files.(tag).in_g6 = in_g6;
        meta.files.(tag).k = k;
        meta.files.(tag).m = m;
        meta.files.(tag).count_in = numel(idx);
        meta.files.(tag).count_classes_in_group = numel(unique(canon));
    end
end


function A = ts2adj_tanner_nopad(H, TS)
% Build adjacency of the induced Tanner subgraph for TS, with ONLY incident CNs.
% Output A is (k+m) x (k+m) with VNs first then CNs.

    if ~islogical(H), H = (H ~= 0); end

    TS = unique(TS(:).');
    subH = H(:, TS);
    subH(sum(subH,2)==0,:) = [];     % keep only CNs that touch TS

    k = numel(TS);
    m = size(subH,1);

    A = false(k+m, k+m);

    [r,c] = find(subH);              % r: CN (1..m), c: VN (1..k)
    idxVN = c;
    idxCN = k + r;

    A(sub2ind([k+m, k+m], idxVN, idxCN)) = true;
    A(sub2ind([k+m, k+m], idxCN, idxVN)) = true;
end


function g6 = adj2graph6(A)
    n = size(A,1);
    if n > 62, error("This simple encoder only supports n <= 62"); end
    
    % Header: n + 63
    header = char(n + 63);
    
    % Get upper triangle bits in order (0,1), (0,2), (1,2), (0,3)...
    % Nauty/Graph6 order: column by column, upper triangle
    bits = [];
    for j = 2:n
        for i = 1:j-1
            bits(end+1) = A(i,j);
        end
    end
    
    % Pad with zeros to a multiple of 6
    padding = 6 - mod(numel(bits), 6);
    if padding < 6
        bits = [bits, zeros(1, padding)];
    end
    
    % Convert to characters
    num_chars = numel(bits) / 6;
    body = setstr(zeros(1, num_chars));
    for k = 1:num_chars
        chunk = bits((k-1)*6 + (1:6));
        % Standard graph6: bit 1 is 2^5, bit 6 is 2^0
        val = sum(chunk .* [32 16 8 4 2 1]);
        body(k) = char(val + 63);
    end
    g6 = string([header body]);
end



function plot_ts_from_H(H, TSrep, opts)
% Non-bipartite trapping-set visualization from H and VN indices.
% Draws:
%  - VNs as circles (main nodes)
%  - CNs as squares near their VN neighborhoods (all CNs shown)
%  - CNs unsatisfied (odd degree within induced subgraph) as filled black squares
%  - CNs satisfied as empty squares
%
% Inputs:
%   H     : parity-check matrix (logical or 0/1)
%   TSrep : vector of VN indices (1-based)
%   opts  : struct with fields:
%       .layout       : 'force' (default), 'circle', 'layered'
%       .show_vn_lbl  : true (default), false
%       .show_cn_lbl  : false (default), true
%       .unsat_def    : 'odd' (default) or 'deg1'
%
% Output: none (creates a figure)

    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'layout'),      opts.layout = "force"; end
    if ~isfield(opts,'show_vn_lbl'), opts.show_vn_lbl = true; end
    if ~isfield(opts,'show_cn_lbl'), opts.show_cn_lbl = false; end
    if ~isfield(opts,'unsat_def'),   opts.unsat_def = "odd"; end

    if ~islogical(H), H = (H ~= 0); end
    TSrep = unique(TSrep(:).');
    k = numel(TSrep);

    % Induced subgraph incidence: keep only CNs that touch TSrep
    subH = H(:, TSrep);
    keep = (sum(subH,2) > 0);
    subH = subH(keep,:);
    m = size(subH,1);

    % VN-CN incidence (k x m)
    B = subH.';   % VNs x CNs

    % CN degrees within induced subgraph
    deg_c = sum(B,1);

    switch lower(opts.unsat_def)
        case 'deg1'
            unsatC = find(deg_c == 1);
        otherwise % 'odd'
            unsatC = find(mod(deg_c,2) == 1);
    end
    satC = setdiff(1:m, unsatC);

        % ---- PLOT: compute layout from VN-CN bipartite graph, then draw manually ----
    figure; ax = axes; hold(ax,'on');

    % Build bipartite adjacency for layout:
    % nodes 1..k are VNs, nodes k+1..k+m are CNs
    A = false(k+m, k+m);
    for c = 1:m
        vs = find(B(:,c)); % 1..k
        A(vs, k+c) = true;
        A(k+c, vs) = true;
    end
    Gb = graph(A);

    % Get coordinates from MATLAB layout engine (do NOT show its edges)
    switch lower(opts.layout)
        case 'circle'
            ptmp = plot(Gb,'Layout','circle','Parent',ax,'Visible','off');
        case 'layered'
            ptmp = plot(Gb,'Layout','layered','Parent',ax,'Visible','off');
        otherwise
            ptmp = plot(Gb,'Layout','force','Parent',ax,'Visible','off');
    end
    Xall = ptmp.XData(:);
    Yall = ptmp.YData(:);
    delete(ptmp);

    Xv = Xall(1:k);          Yv = Yall(1:k);
    Xc0 = Xall(k+1:k+m);     Yc0 = Yall(k+1:k+m);

    % Optional: tiny deterministic jitter for CNs to reduce overlap (keep small)
    for c = 1:m
        ang = 2*pi*mod(c*0.61803398875,1);
        r   = 0.04; % much smaller than before since layout already separates
        Xc0(c) = Xc0(c) + r*cos(ang);
        Yc0(c) = Yc0(c) + r*sin(ang);
    end

      % --- polish: center + PCA-rotate for a more "symmetric" look ---
    XY = [Xall(:), Yall(:)];
    XY = XY - mean(XY,1);                % center

    % PCA rotation (align principal axis to x)
    C = cov(XY);
    [V,D] = eig(C);
    [~,ix] = max(diag(D));
    v1 = V(:,ix);
    ang = atan2(v1(2), v1(1));
    R = [cos(-ang) -sin(-ang); sin(-ang) cos(-ang)];
    XY = (R*XY.').';

    Xall = XY(:,1); Yall = XY(:,2);

    Xv  = Xall(1:k);        Yv  = Yall(1:k);
    Xc0 = Xall(k+1:k+m);    Yc0 = Yall(k+1:k+m);

    % small deterministic jitter for CNs (keep subtle)
    for c = 1:m
        angj = 2*pi*mod(c*0.61803398875,1);
        rj   = 0.03;
        Xc0(c) = Xc0(c) + rj*cos(angj);
        Yc0(c) = Yc0(c) + rj*sin(angj);
    end

    % --- compute label offsets in data units ---
    xmin = min(Xall); xmax = max(Xall);
    ymin = min(Yall); ymax = max(Yall);
    sx = max(xmax - xmin, eps);
    sy = max(ymax - ymin, eps);
    dx = 0.02 * sx;   % label offset
    dy = 0.02 * sy;

    % ---- draw VNs (white circles, black stroke) ----
    vnSize = 95;
    scatter(ax, Xv, Yv, vnSize, 'o', ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineWidth', 1.5);

    if opts.show_vn_lbl
        for i = 1:k
            text(ax, Xv(i)+dx, Yv(i)+dy, sprintf('%d', TSrep(i)), ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','bottom', ...
                'FontSize', 9, ...
                'Color', [0 0 0], ...
                'Interpreter','none');
        end
    end

    % ---- draw CNs + spokes (elegant styling) ----
    cnSize = 175;
    spokeColor = [0.72 0.72 0.72];

    for c = 1:m
        cx2 = Xc0(c); cy2 = Yc0(c);
        vs = find(B(:,c));

        if ismember(c, unsatC)
            face = [0 0 0]; edge = [0 0 0]; lw = 1.8;
        else
            face = [1 1 1]; edge = [0 0 0]; lw = 1.2;
        end

        scatter(ax, cx2, cy2, cnSize, 's', ...
            'MarkerFaceColor', face, ...
            'MarkerEdgeColor', edge, ...
            'LineWidth', lw);

        % spokes
        for v = vs(:).'
            plot(ax, [cx2 Xv(v)], [cy2 Yv(v)], '-', ...
                'Color', spokeColor, 'LineWidth', 1.0);
        end

        if opts.show_cn_lbl
            text(ax, cx2+dx, cy2+dy, sprintf('c%d', c), ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','bottom', ...
                'FontSize', 9, ...
                'Color', [0 0 0], ...
                'Interpreter','none');
        end
    end

    axis(ax,'equal'); axis(ax,'off');

end

function [isIso, key1, key2, info] = are_iso_ts_labelg(H, TS1, TS2, nautyDirLinux)
    if ~islogical(H), H = (H ~= 0); end
    TS1 = unique(TS1(:).'); TS2 = unique(TS2(:).');

    % induced CN counts
    subH1 = H(:,TS1); subH1(sum(subH1,2)==0,:) = [];
    subH2 = H(:,TS2); subH2(sum(subH2,2)==0,:) = [];
    k1 = numel(TS1);  m1 = size(subH1,1);
    k2 = numel(TS2);  m2 = size(subH2,1);

    info = struct('k1',k1,'m1',m1,'k2',k2,'m2',m2);

    % Must match (k,m) for colored isomorphism
    if k1 ~= k2 || m1 ~= m2
        isIso = false; key1 = ""; key2 = "";
        return;
    end

    labelg = string(nautyDirLinux) + "/labelg";

    outDir = "C:\temp\nauty_tmp";
    if ~isfolder(outDir), mkdir(outDir); end
    in_g6 = fullfile(outDir, sprintf("pair_k%d_m%d.g6", k1, m1));

    fid = fopen(in_g6,"w");
    if fid<0, error("Cannot open %s", in_g6); end

    A1 = ts2adj_tanner_nopad(H, TS1);
    A2 = ts2adj_tanner_nopad(H, TS2);
    fprintf(fid, "%s\n", adj2graph6(A1));
    fprintf(fid, "%s\n", adj2graph6(A2));
    fclose(fid);

    fopt = sprintf('-fa^%db^%d', k1, m1);
    cmd = sprintf('wsl bash -lc "%s -q %s \\"%s\\" 2>/dev/null"', ...
        labelg, fopt, win2wsl_path(in_g6));

    [status,out] = system(cmd);
    if status~=0
        error("labelg failed:\n%s\nCMD:\n%s", out, cmd);
    end

    canon = string(splitlines(string(out)));
    canon = strip(canon);
    canon = canon(canon ~= "");
    if numel(canon) ~= 2
        error("Expected 2 output lines from labelg, got %d", numel(canon));
    end

    key1 = canon(1); key2 = canon(2);
    isIso = (key1 == key2);
end

function inv = ts_invariants(H, TS)
    if ~islogical(H), H = (H ~= 0); end
    TS = unique(TS(:).');
    subH = H(:,TS);
    subH(sum(subH,2)==0,:) = [];
    B = subH.';                 % k x m incidence
    dv = sort(sum(B,2));        % VN degrees within induced
    dc = sort(sum(B,1)).';      % CN degrees within induced
    inv = struct('k',numel(TS),'m',size(subH,1), ...
                 'dv',dv(:).', 'dc',dc(:).');
end

function A = graph6_to_adj(g6)
    s = char(g6);
    n = double(s(1)) - 63;
    A = false(n, n);
    
    data = double(s(2:end)) - 63;
    bits = [];
    for i = 1:numel(data)
        chunk_bits = bitget(data(i), 6:-1:1);
        bits = [bits, chunk_bits];
    end
    
    % Fill upper triangle in column-major order
    idx = 1;
    for j = 2:n
        for i = 1:j-1
            if bits(idx)
                A(i,j) = true;
                A(j,i) = true;
            end
            idx = idx + 1;
        end
    end
end

function show_incidence(H, TSrep)
    if ~islogical(H), H = (H ~= 0); end
    TSrep = unique(TSrep(:).');

    subH = H(:, TSrep);
    keep = sum(subH,2) > 0;
    subH = subH(keep,:);   % m x k

    figure;
    imagesc(subH);
    colormap(gray);
    axis tight;
    xlabel('VN (order in TSrep)');
    ylabel('CN (local index)');
    title(sprintf('Incidence m=%d, k=%d', size(subH,1), size(subH,2)));
end