function Lp = enumeratePaths(NL,CL,amax,K,g)

Lp = cell(size(NL,1),amax);
for StartingNode = 1 : size(NL,2);
    Lp{StartingNode}{1}{1} = NL{StartingNode};
    Lp{StartingNode}{1}{2} = CL{StartingNode};
    for k=2 : min(K,amax) 
        
        if mod((2*k)/2-1,2)==1
            [lista_percorsiV,lista_percorsiC]=FindPaths_Bani(StartingNode,k/2,NL,CL,0); 
        else
            [lista_percorsiV,lista_percorsiC]=FindPaths_Bani(StartingNode,(k-1)/2,NL,CL,1);
        end
            Lp{StartingNode}{k}{1} = lista_percorsiV;
            Lp{StartingNode}{k}{2} = lista_percorsiC;
    end
end
end
        