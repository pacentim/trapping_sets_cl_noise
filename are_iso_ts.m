function [isIso, key1, key2] = are_iso_ts(H, TS1, TS2, nautyDirLinux)
    if ~islogical(H), H = (H ~= 0); end

    TS1 = unique(TS1(:).');  TS2 = unique(TS2(:).');
    k1 = numel(TS1);         k2 = numel(TS2);
    if k1 ~= k2
        isIso = false; key1=""; key2=""; return;
    end
    k = k1;

    % induced CN counts
    subH1 = H(:,TS1); subH1 = subH1(any(subH1,2),:); m1 = size(subH1,1);
    subH2 = H(:,TS2); subH2 = subH2(any(subH2,2),:); m2 = size(subH2,1);
    Mk = max(m1,m2);

    % build padded bipartite adjacencies (n = k + Mk for both)
    A1 = ts2adj_tanner_padded(H, TS1, Mk);
    A2 = ts2adj_tanner_padded(H, TS2, Mk);

    g61 = adj2graph6(A1);
    g62 = adj2graph6(A2);

    % write temporary file with two graphs
    in_g6 = fullfile(tempdir, "two_ts.g6");
    fid = fopen(in_g6,"w");
    fprintf(fid, "%s\n%s\n", g61, g62);
    fclose(fid);

    labelg = nautyDirLinux + "/labelg";
    fopt = sprintf('-fa^%db^%d', k, Mk);

    cmd = sprintf('wsl bash -lc "%s -q %s %s"', labelg, fopt, win2wsl_path(in_g6));
    [status,out] = system(cmd);
    if status ~= 0
        error("labelg failed:\n%s", out);
    end

    canon = string(splitlines(string(out)));
    canon = canon(canon ~= "");
    key1 = canon(1);
    key2 = canon(2);

    isIso = (key1 == key2);
end

function A = ts2adj_tanner_padded(H, TS, Mk)
    TS = unique(TS(:).');
    subH = H(:, TS);
    subH = subH(any(subH,2),:);

    k = numel(TS);
    m = size(subH,1);
    if m > Mk, error("m > Mk"); end

    A = false(k+Mk, k+Mk);
    [r,c] = find(subH);      % r CN in 1..m, c VN in 1..k
    A(sub2ind([k+Mk,k+Mk], c, k+r)) = true;
    A(sub2ind([k+Mk,k+Mk], k+r, c)) = true;
end
