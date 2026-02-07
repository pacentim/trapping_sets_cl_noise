function hFig = draw_ts_induced_graph(Hsub, varargin)
%DRAW_TS_INDUCED_GRAPH Draw Tanner subgraph induced by a trapping set.
%
%   hFig = draw_ts_induced_graph(Hsub, 'Name', Value, ...)
%
% Input
%   Hsub : (m x n) binary submatrix for the induced subgraph:
%          rows = check nodes (CNs), cols = variable nodes (VNs).
%          Entry 1 means an edge between CN_i and VN_j.
%
% Plot style
%   Normal (non-bipartite) layout using a force-directed embedding.
%   VNs are circles, CNs are squares. Unsatisfied CNs (odd degree) are highlighted.
%
% Name-Value options
%   'ShowLabels'   : true/false (default true)
%   'VNLabels'     : string/cellstr length n (default "v1"...)
%   'CNLabels'     : string/cellstr length m (default "c1"...)
%   'Layout'       : 'force'|'layered'|'circle' etc. (default 'force')
%   'FigureName'   : figure name (default 'Trapping set subgraph')
%   'NodeSize'     : marker size (default 80)
%   'EdgeAlpha'    : edge transparency (default 0.25)
%   'UnsatCNColor' : RGB (default [0.85 0.25 0.25])
%
% Output
%   hFig : figure handle

    p = inputParser;
    p.addParameter('ShowLabels', true, @(x)islogical(x) || isnumeric(x));
    p.addParameter('VNLabels', [], @(x) isstring(x) || iscellstr(x) || isempty(x));
    p.addParameter('CNLabels', [], @(x) isstring(x) || iscellstr(x) || isempty(x));
    p.addParameter('Layout', 'force', @(s)ischar(s) || isstring(s));
    p.addParameter('FigureName', 'Trapping set subgraph', @(s)ischar(s) || isstring(s));
    p.addParameter('NodeSize', 80, @(x)isnumeric(x) && isscalar(x) && x>0);
    p.addParameter('EdgeAlpha', 0.25, @(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    p.addParameter('UnsatCNColor', [0.85 0.25 0.25], @(x)isnumeric(x) && numel(x)==3);
    p.parse(varargin{:});
    opt = p.Results;

    Hsub = logical(Hsub);
    [mC, nV] = size(Hsub);

    % Unsatisfied CNs: odd degree within induced subgraph
    degC = full(sum(Hsub, 2));
    unsatCN = mod(degC, 2) == 1;

    % Build a plain (non-bipartite) undirected graph with node set [VNs, CNs]
    % Node ids: 1..nV are VNs, nV+1..nV+mC are CNs
    [ci, vj] = find(Hsub);      % ci in 1..mC, vj in 1..nV
    s = vj(:);
    t = nV + ci(:);
    G = graph(s, t);

    % Labels
    if isempty(opt.VNLabels)
        vnNames = "v" + string(1:nV);
    else
        vnNames = string(opt.VNLabels(:))';
        if numel(vnNames) ~= nV, error('VNLabels must have length nV.'); end
    end
    if isempty(opt.CNLabels)
        cnNames = "c" + string(1:mC);
    else
        cnNames = string(opt.CNLabels(:))';
        if numel(cnNames) ~= mC, error('CNLabels must have length mC.'); end
    end
    nodeNames = [vnNames, cnNames];
    G.Nodes.Name = nodeNames(:);

    % Figure
    hFig = figure('Name', char(opt.FigureName), 'Color', 'w');
    ax = axes('Parent', hFig); %#ok<LAXES>
    hold(ax, 'on');
    axis(ax, 'off');

    % Plot edges only, get coordinates from layout
    ph = plot(ax, G, ...
        'Layout', char(opt.Layout), ...
        'NodeLabel', {}, ...
        'Marker', 'none');
    ph.EdgeAlpha = opt.EdgeAlpha;

    X = ph.XData(:);
    Y = ph.YData(:);

    % Scatter nodes with different shapes
    vnIdx = 1:nV;
    cnIdx = (nV+1):(nV+mC);

    % VN: circles
    scatter(ax, X(vnIdx), Y(vnIdx), opt.NodeSize, ...
        'o', 'filled', ...
        'MarkerFaceAlpha', 1, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5);

    % CN: squares (satisfied vs unsatisfied)
    cnSat = cnIdx(~unsatCN(:)');
    cnUns = cnIdx(unsatCN(:)');

    if ~isempty(cnSat)
        scatter(ax, X(cnSat), Y(cnSat), opt.NodeSize, ...
            's', 'filled', ...
            'MarkerFaceAlpha', 1, ...
            'MarkerFaceColor', [0.25 0.25 0.25], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.5);
    end
    if ~isempty(cnUns)
        scatter(ax, X(cnUns), Y(cnUns), opt.NodeSize, ...
            's', 'filled', ...
            'MarkerFaceAlpha', 1, ...
            'MarkerFaceColor', opt.UnsatCNColor, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8);
    end

    % Optional labels
    if opt.ShowLabels
        % Slight offset so labels don't sit on markers (portable replacement for range)
        rx = max(X) - min(X);
        ry = max(Y) - min(Y);
        if rx == 0, rx = 1; end
        if ry == 0, ry = 1; end
        dx = 0.01 * rx;
        dy = 0.01 * ry;

        for u = 1:(nV+mC)
            text(ax, X(u)+dx, Y(u)+dy, nodeNames(u), ...
                'Interpreter', 'none', ...
                'FontSize', 9);
        end
    end

    title(ax, sprintf('|V|=%d, |C|=%d, unsat CNs=%d, edges=%d', ...
        nV, mC, nnz(unsatCN), numedges(G)), 'Interpreter', 'none');

    hold(ax, 'off');
end
