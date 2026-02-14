outDir = "ts_plots";
if ~isfolder(outDir)
    mkdir(outDir);
end

load('failure_sets.mat');
nShow = length(TS);

for i = 1:nShow
    
TSrep = TS{i};

    % Create invisible figure
    fig = figure('Visible','off');
    
    plot_ts_from_H(H, TSrep, struct( ...
        'show_vn_lbl', true, ...
        'show_cn_lbl', true, ...
        'unsat_def', "odd"));
    
    % Build filename
    fname = fullfile(outDir, sprintf("TS_%03d.png", i));
    
    % Save at high resolution
    exportgraphics(fig, fname, 'Resolution', 300);
    
    % Close to avoid memory accumulation
    close(fig);
end