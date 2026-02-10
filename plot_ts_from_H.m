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
