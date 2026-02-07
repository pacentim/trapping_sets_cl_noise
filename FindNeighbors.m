function [NL, CL] = FindNeighbors(H)

NL = cell(1,size(H,2));
CL = cell(1,size(H,2));

for i=1:size(H,2) %per tutte le colonne di H
    Rows = find(H(:, i)); % indici di ogni colonna
    NeighborList = NaN*ones(1e3,1);
    CheckNodeList = NaN*ones(1e3,1);
    count = 1;
    for Row = Rows' %Row è un vettore riga che contiene gli indici della iesima colonna
        newNeighbors = find(H(Row, :))'; %i neighbors di ogni nodo check
        CheckNodeList(count:count+length(newNeighbors)-1,1) = ones(length(newNeighbors),1) * Row;  %x ogni neighbor è indicato l'indice del nodo check corrispondente
        NeighborList(count:count+length(newNeighbors)-1,1) = newNeighbors;
        count = count + length(newNeighbors);
    end
    CheckNodeList = CheckNodeList(~isnan(CheckNodeList));
    NeighborList = NeighborList(~isnan(NeighborList));
    %%%%%%%%%%%%
    SameNodeIdx = NeighborList == i; %elimino i falsi neighbors che corrispondono al nodo check stesso
    %%%%%%%%%%%%%
    %SameNodeIdx = find(NeighborList == VarNodeIdx);
    NeighborList(SameNodeIdx) = [];
    CheckNodeList(SameNodeIdx) = [];
    NL{i}=NeighborList;
    CL{i}=CheckNodeList;
end


