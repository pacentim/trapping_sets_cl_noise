function [ PRIMES, finalPOSjustV] = EnumerateCycles_Bani_NON_QC(H, CyclesLength,Lp)

close all;
POS={};         %Matrice Generale contenente le posizioni degli 1 dei loop
finalPOS=cell(7,1);    %Matrice finale contentente le posizioni degli 1 solamente dei loop che non si ripetono
PRIMES = {};
 finalPOSjustV = [];
ind=1;          %Indice Matrice vettori posizioni_loop finale
ind2=1;
%ind_clone=1;    %Indice Matrice cloni trovati

checking_mat=zeros(0,2*CyclesLength);

[n,m]=size(H);  %Acquisizione delle dimensioni della matrice H

T=1:m;

% Rappresentazione della matrice H.

[row_idxs,col_idxs]=find(H);
coord=[row_idxs,col_idxs];

% Selezione delle colonne di interesse
C=T;
H1=H(:,C);
cont_fPOS=1;
%spy(H1);
%Ciclo Principale
for KK=1:size(H1,2)
    %for KK=length(C):-1:1
    StartingNode = (KK);  % nodo di partenza
    
    % Valutazione dell'ordine in cui sono toccate le colonne
        lista_percorsiV = Lp{StartingNode}{CyclesLength/2}{1};
        lista_percorsiC = Lp{StartingNode}{CyclesLength/2}{2};
    
    %se la lunghezza del ciclo è 2*k, non mi serve prolungare al check
    if mod(CyclesLength/2-1,2)==1

        %[lista_percorsiV,lista_percorsiC]=FindPaths_Bani(StartingNode,CyclesLength/4,NL,CL,0);
        for i=1:size(lista_percorsiV,1)
            LoopsIdxs = lista_percorsiV(:, end) == lista_percorsiV(i,end);
            LoopsIdxs(1:i)=0;
            Tlista_percorsiV = lista_percorsiV(LoopsIdxs, :);
            Tlista_percorsiC = lista_percorsiC(LoopsIdxs, :);
            
            p1V=lista_percorsiV(i,2:(end-1));
            p1C=lista_percorsiC(i,2:(end));
            
            for J=1:size(Tlista_percorsiV,1)
                Tp1V=Tlista_percorsiV(J,2:(end-1));
                Tp1C=Tlista_percorsiC(J,2:(end));
                if isempty(intersect(p1V,Tp1V)) && isempty(intersect(Tp1C,p1C))
                    candidate=union(lista_percorsiV(i,:),Tlista_percorsiV(J,:));
                    finalPOSjustV{1,ind}=sort(candidate);
                    ind=ind+1;
                    
                end
            end
        end
    %se la lunghezza del ciclo è 2*k, mi serve prolungare al check
    else
%         [lista_percorsiV,lista_percorsiC]=FindPaths_Bani(StartingNode,(CyclesLength-2)/4,NL,CL,1);
        
        for i=1:size(lista_percorsiV,1)
            LoopsIdxs = lista_percorsiC(:, end) == lista_percorsiC(i,end);
            LoopsIdxs(1:i)=0;
            Tlista_percorsiV = lista_percorsiV(LoopsIdxs, :);
            Tlista_percorsiC = lista_percorsiC(LoopsIdxs, :);
            
            p1V=lista_percorsiV(i,2:end);
            p1C=lista_percorsiC(i,1:(end-1));
            
            for J=1:size(Tlista_percorsiV,1)
                Tp1V=Tlista_percorsiV(J,2:end);
                Tp1C=Tlista_percorsiC(J,1:(end-1));
                if isempty(intersect(p1V,Tp1V)) && isempty(intersect(Tp1C,p1C))
                    candidate=union(lista_percorsiV(i,:),Tlista_percorsiV(J,:));
                        H_par=H(:,candidate);
                        somma_Hpar=sum(H_par,2);
                        unsati=sort(find(mod(somma_Hpar,2)==1 & somma_Hpar>0)');
                        b=length(unsati);
                        finalPOSjustV{1,ind}=sort(candidate);
                        ind=ind+1;
                    end
                end
            end
            
        end
    end
    
    
end

