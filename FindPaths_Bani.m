function [PathList,CheckNodePathList] = FindPaths_Bani( VarNodeIdx, PathDepth,NL,CL,EndCheck)
%%%%%%%%%%


PathList = VarNodeIdx;
CheckNodePathList = -1;
for currDepth = 1:PathDepth
    %tempPathList = [];
    %tempCheckNodePathList = [];
    PathNr = size(PathList, 1);
    %%%%%%%%%%%%
%     if PathDepth<6
        tempPathList=zeros(3*5^(PathDepth)*2^(PathDepth-1)+1,currDepth+1);
        tempCheckNodePathList=tempPathList;
%     else
%         tempPathList = [];
%         tempCheckNodePathList = [];
%     end
    cont=1;
    %%%%%%%%%%%%%%%
    for I = 1:PathNr
        currNode = PathList(I, end);
        
        %%%%%%%%%%%
        NeighborList=NL{currNode};
        CheckNodeList=CL{currNode};
        %%%%%%%%%%%%%
        
        newPathList = [ones(size(NeighborList,1), 1)*PathList(I,:) NeighborList];
        newCheckNodePathList = [ones(size(CheckNodeList,1), 1)*CheckNodePathList(I,:) CheckNodeList];
        
        for J=2:currDepth
            id= newCheckNodePathList(:,end)-newCheckNodePathList(:,J)==0;
            newPathList(id,:) = [];
            newCheckNodePathList(id,:) = [];
        end
        
%         if PathDepth<6
            snp=size(newPathList,1);
            tempPathList(cont:cont+snp-1,:)=newPathList;
            tempCheckNodePathList(cont:cont+snp-1,:)=newCheckNodePathList;
            cont=cont+snp;
%         else
%             tempPathList = [tempPathList; newPathList];
%             tempCheckNodePathList = [tempCheckNodePathList; newCheckNodePathList];
%         end
    end
    
    %%%%%%%%%%%%%%%%%%%
    ss=~any(tempPathList,2);
    tempPathList( ss, : ) = [];
    tempCheckNodePathList( ss, : ) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%
   
    PathList = tempPathList;
    CheckNodePathList = tempCheckNodePathList;
    
    if EndCheck && currDepth==PathDepth
        cont=1;
        for I = 1:size(PathList,1)
            currNode = PathList(I, end);
            
            %%%%%%%%%%%
            %NeighborList=NL{currNode};
            CheckNodeList=unique(CL{currNode},'rows');
            %%%%%%%%%%%%%
            
            newPathList = [ones(size(CheckNodeList,1), 1)*PathList(I,:)];
            newCheckNodePathList = [ones(size(CheckNodeList,1), 1)*CheckNodePathList(I,:) CheckNodeList];
            
            for J=2:currDepth+1
                id= newCheckNodePathList(:,end)-newCheckNodePathList(:,J)==0;
                newPathList(id,:) = [];
                newCheckNodePathList(id,:) = [];
            end
            
%             if PathDepth<6
                snp=size(newPathList,1);
                tempPathList(cont:cont+snp-1,:)=newPathList;
                tempCheckNodePathList(cont:cont+snp-1,:)=newCheckNodePathList(:,2:end);
                cont=cont+snp;
%             else
%                 tempPathList = [tempPathList; newPathList];
%                 tempCheckNodePathList = [tempCheckNodePathList; newCheckNodePathList];
%             end
        end
        %%%%%%%%%%%%%%%%%%%
        ss=~any(tempPathList,2);
        tempPathList( ss, : ) = [];
        tempCheckNodePathList( ss, : ) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        PathList = tempPathList;
        CheckNodePathList = tempCheckNodePathList;
    end
    
end

