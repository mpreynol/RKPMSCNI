classdef Cloud < handle
    %CLOUD is a list of node objects for an RKPM solution space
    % Currently using 2D
    
    properties
        nodalData=[]; % Array of nodal Cordinates
        Nodes=[]; % List of Object Nodes
        numberOfNodes; % Number of nodes in the cloud
        order; % Order of basis Functions
        aArray; % Array of a values
        aScale; % General Scale for the dilation parameter
        u; % Solutions *** Lack of Kronecker Means this isn't the PDE Solution Yet
        Verts; % Voroni Celll Vertex
        Cells; % Cell Structure of the Voroni Diagram
    end
    
    methods
        function obj=Cloud(nodalData,order,aScale,Verts,Cells)
            %Constructor
            obj.nodalData=nodalData;
            obj.numberOfNodes=size(nodalData,1);
            obj.order=order;
            obj.aScale=aScale;
            obj.Verts=Verts;
            obj.Cells=Cells;
            obj.aArray=obj.findDilation();
            obj.aArray=obj.aArray*aScale;
            obj.setNodeObjects();
        end
        
        function setNodeObjects(obj)
            obj.Nodes=RKNode.empty(obj.numberOfNodes,0);
            for i=1:obj.numberOfNodes
                VoroniCords=obj.Verts(obj.Cells{i}(:),:);
                obj.Nodes(i)=RKNode(obj.nodalData(i,1),[obj.nodalData(i,2);obj.nodalData(i,3)],obj.aArray(i),obj.order,[obj.nodalData(i,4);obj.nodalData(i,5)],VoroniCords,obj);
            end
        end
        
        function A=distToAll(obj)
            % Matrix with each element a data structure with distance vector between points i and j (sym).
            A=zeros(obj.numberOfNodes,obj.numberOfNodes,3);
            for i=1:obj.numberOfNodes
                for j=1:obj.numberOfNodes
                    if i~=j
                        A(i,j,1)=obj.nodalData(i,2)-obj.nodalData(j,2);
                        A(i,j,2)=obj.nodalData(i,3)-obj.nodalData(j,3);
                        A(i,j,3)=sqrt(A(i,j,1)^2+A(i,j,2)^2);
                    end
                end
            end
        end
        
        function aArray=findDilation(obj)
            % Method returns an optimal 'a' value for the mesh: smallest
            % nonCoplaner Point
            A=obj.distToAll();
            aArray=ones(obj.numberOfNodes,1)*inf;
            for i=1:obj.numberOfNodes
               for j=1:obj.numberOfNodes
                   if i~=j
                      if A(i,j,1)~=0 && A(i,j,2)~=0 % Point is not CoPlanar
                         if A(i,j,3)<aArray(i)
                             aArray(i)=A(i,j,3);
                         end
                      end
                   end
               end
            end
        end
        
        function plotCloud(obj)
            for i=1:obj.numberOfNodes
                refresh
                if obj.nodalData(i,5)==0
                    scatter(obj.nodalData(i,2),obj.nodalData(i,3),'ro','filled')   
                else
                    scatter(obj.nodalData(i,2),obj.nodalData(i,3),'ko','filled')
                end
                hold on
            end
        end
        
        function plotVoroni(obj,scale)
            %Method Plots a Voroni Diagram
            % Scale is the Scale For Plotting the Normal
            for i=1:obj.numberOfNodes
                localNode=obj.Nodes(i);
                Points=localNode.Voron.VoroniCords;
                Points=[Points;Points(1,:)];
                for j=1:size(Points,1)-1
                    plot([Points(j,1) Points(j+1,1)], [Points(j,2) Points(j+1,2)],'b');
                    hold on
                end
                for j=1:size(localNode.Voron.midPoint,1)
                   %plot(localNode.midPoint(j,1),localNode.midPoint(j,2),'r*')
                   r=[localNode.Voron.normalPoint(j,1)-localNode.Voron.midPoint(j,1), ...
                      localNode.Voron.normalPoint(j,2)-localNode.Voron.midPoint(j,2)];
                  r=r*scale;
                  plot([localNode.Voron.midPoint(j,1) localNode.Voron.midPoint(j,1)+r(1)], ...
                      [localNode.Voron.midPoint(j,2) localNode.Voron.midPoint(j,2)+r(2)],'r-')
                  %pause(0.1)
                end
            end
        end
        
        function plotVoroniDeformed(obj,scale)
            %Method Plots a Voroni Diagram in deformed condition
            % Scale is the Scale For Plotting the Normal
            for i=1:obj.numberOfNodes
                localNode=obj.Nodes(i);
                Points=localNode.Voron.VoroniCords;
                for j=1:size(localNode.Voron.normal,1)
                    disp=obj.returnInterpolatedU(localNode.Voron.VoroniCords(j,:)');
                    Points(j,:)=Points(j,:)+scale*disp';
                end  
                Points=[Points;Points(1,:)];
                for j=1:size(Points,1)-1
                    plot([Points(j,1) Points(j+1,1)], [Points(j,2) Points(j+1,2)],'-o');
                    hold on
                end
            end
        end
        
        
        function [ufull] = reAssembleUnknowns(obj,ureduced,BE)
            % Method reassembles a full 'u' vector for a reduced oned
            L=BE==-inf;
            ufull=zeros(length(BE),1);
            counter=1;
            for i=1:length(BE)
                if L(i)==1
                    ufull(i)=ureduced(counter);
                    counter=counter+1;
                else
                    ufull(i)=BE(i);
                end
            end
        end
        
        function parseSolution(obj,u)
            % Method insets solved for coeff's 'u' into the Cloud and Nodes
            obj.u=u;
            for i=1:obj.numberOfNodes
                dof=[obj.Nodes(i).nodeNumber*2-1;obj.Nodes(i).nodeNumber*2];
                obj.Nodes(i).u=u(dof);
            end
        end
        
        function d=plotCloudU(obj,scale)
            d=zeros(2,obj.numberOfNodes);
            for j=1:obj.numberOfNodes
                Cords=obj.Nodes(j).cordinates;
                deformation=zeros(2,1);
                for i=1:obj.numberOfNodes
                    if Cords(1)-obj.Nodes(i).cordinates(1)<=obj.Nodes(i).a && Cords(2)-obj.Nodes(i).cordinates(2)<=obj.Nodes(i).a
                        deformation=deformation+obj.Nodes(i).sF.getValue([Cords(1);Cords(2)])*obj.Nodes(i).u;
                    end
                end
                d(:,j)=deformation;
                scatter(Cords(1)+scale*deformation(1),Cords(2)+scale*deformation(2),'bo')
                hold on
            end
        end
        
        function [K,F]=integrateDomain(obj,C,Q)
            n=obj.numberOfNodes;
            K=zeros(n*2,n*2);
            F=zeros(n*2,1);
            for k=1:n % Loop through Nodes
                obj.Nodes(k).CellDeriv.setData(Q);
                k/n*100
                A=obj.Nodes(k).Voron.area;
                for a=1:n % Loop over Shape Functions (Rows)
                    if obj.Nodes(k).CellDeriv.B(a,1)~=0 % 
                        b_a1=obj.Nodes(k).CellDeriv.B(a,1); % Get a's Component of the B Matrix
                        b_a2=obj.Nodes(k).CellDeriv.B(a,2); % Get a's Component of the B Matrix
                        Ba=[b_a1,0;0,b_a2;b_a2,b_a1];
                        for b=1:n % Loop over Shape Functions (Columns)
                            if b>=a % Only assemble the upper half of matrix
                                if obj.Nodes(k).CellDeriv.B(b,1)~=0 % 
                                    b_b1=obj.Nodes(k).CellDeriv.B(b,1); % Get b's Component of the B Matrix
                                    b_b2=obj.Nodes(k).CellDeriv.B(b,2); % Get b's Component of the B Matrix
                                    Bb=[b_b1,0;0,b_b2;b_b2,b_b1];
                                    Kab=Ba'*C*Bb*A;
                                    K(2*a-1,2*b-1)=K(2*a-1,2*b-1)+Kab(1,1);
                                    K(2*a-1,2*b)=K(2*a-1,2*b)+Kab(1,2);
                                    K(2*a,2*b-1)=K(2*a,2*b-1)+Kab(2,1);
                                    K(2*a,2*b)=K(2*a,2*b)+Kab(2,2);
                                end
                            end
                        end
                    end
                end
                F=obj.Nodes(k).CellDeriv.F+F;
            end
            K=triu(K)+tril(K',-1);
        end
        
        function F=integrateBoundary(obj,Mesh,BN)
            F=zeros(2*obj.numberOfNodes,1);
            for j=1:Mesh.noElements
                h=BN(Mesh.Elements(j).dof);  
                if sum((sum(h~=0))) % Then we have tractions on the element
                    nInt=Mesh.Elements(j).orderInt;
                    G1=Mesh.Elements(j).G1;
                    if (h(1)~=0 && h(3)~=0) || (h(2)~=0 && h(4)~=0)
                        % Surface 1: eta=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),-1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            hInt=[h(1);h(2)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(3);h(4)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(3)~=0 && h(5)~=0) || (h(4)~=0 && h(6)~=0)
                        % Surface 2: xi=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            hInt=[h(3);h(4)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(5);h(6)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(5)~=0 && h(7)~=0) || (h(6)~=0 && h(8)~=0)
                        % Surface 3: eta=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            hInt=[h(5);h(6)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(7);h(8)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(1)~=0 && h(7)~=0) || (h(2)~=0 && h(8)~=0)
                        % Surface 4: xi=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(-1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            hInt=[h(7);h(8)]*(1-(1+G1(i,1))/2)+(1+G1(i,1))/2*[h(1);h(2)];
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*hInt*js*G1(i,2);
                                end
                            end
                        end
                    end
                end
            end
            close(hw);
        end
        function F=integrateExactBoundary(obj,Mesh,BN,exactTraction)
            F=zeros(2*obj.numberOfNodes,1);
            hw=waitbar(0/Mesh.noElements,'Computing Boundary Integration');
            for j=1:Mesh.noElements
                h=BN(Mesh.Elements(j).dof);
                waitbar(j/Mesh.noElements,hw,'Computing Boundary Integration')
                if sum((sum(h~=0))) % Then we have tractions on the element
                    nInt=Mesh.Elements(j).orderInt;
                    G1=Mesh.Elements(j).G1;
                    if (h(1)~=0 && h(3)~=0) || (h(2)~=0 && h(4)~=0)
                        % Surface 1: eta=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),-1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*exactTraction(Cords(1))*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(3)~=0 && h(5)~=0) || (h(4)~=0 && h(6)~=0)
                        % Surface 2: xi=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*exactTraction(Cords(2))*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(5)~=0 && h(7)~=0) || (h(6)~=0 && h(8)~=0)
                        % Surface 3: eta=1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(G1(i,1),1);
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xxi)^2+(Mesh.Elements(j).Shape.Yxi)^2);
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*exactTraction(Cords(1))*js*G1(i,2);
                                end
                            end
                        end
                    end
                    if (h(1)~=0 && h(7)~=0) || (h(2)~=0 && h(8)~=0)
                        % Surface 4: xi=-1
                        for i=1:nInt % perform Guass Integration [-1,1] over domain
                            Mesh.Elements(j).Shape.setAll(-1,G1(i,1));
                            Cords=[Mesh.Elements(j).Shape.X;Mesh.Elements(j).Shape.Y];
                            js=sqrt((Mesh.Elements(j).Shape.Xeta)^2+(Mesh.Elements(j).Shape.Yeta)^2);
                            for w=1:obj.numberOfNodes
                                if norm(Cords-obj.Nodes(w).cordinates)<=obj.Nodes(w).a
                                    F(2*w-1:2*w)=F(2*w-1:2*w)+obj.Nodes(w).sF.getValue(Cords)*exactTraction(Cords(2))*js*G1(i,2);
                                end
                            end
                        end
                    end
                end
            end
            close(hw);
        end
        
        function U=returnInterpolatedU(obj,x)
            U=zeros(2,1);
            for i=1:obj.numberOfNodes
                U=U+obj.Nodes(i).sF.getValue(x)*obj.Nodes(i).u;
            end
        end
        
        function test=checkPU(obj)
            % Test for Partial Unity:
            test=0;
            for j=1:obj.numberOfNodes % Loops through nodes
                for g=1:size(obj.Nodes(j).Voron.VoroniCords,1) % Loops through Edges
                    xTest=obj.Nodes(j).Voron.midPoint(g,:);
                    pu=0;
                    for i=1:obj.numberOfNodes
                        pu=pu+obj.Nodes(i).sF.getValue([xTest(1);xTest(2)]);
                    end
                        pu=abs(pu-1);
                    if pu>test
                        test=pu;
                    end
                end
                
            end
        end
        
        function test=checkNU(obj,Mesh)
            % Test for Partial Nullity:
            test=zeros(2,1);
            hw=waitbar(0/Mesh.noElements,'Checking Partial Nullity at Int Points');
            for j=1:Mesh.noElements
                waitbar(j/Mesh.noElements,hw,'Checking Partial Nullity at Int Points');
                for g=1:Mesh.orderInt
                    xTest=Mesh.Elements(j).getIntCord(g);
                    nu=zeros(2,1);
                    for i=1:obj.numberOfNodes
                        nu=nu+obj.Nodes(i).sF.getValueDx([xTest(1);xTest(2)]);
                    end
                    if abs(nu)>test
                        test=abs(nu);
                    end
                end
                
            end
            close(hw)
        end
        
        function test=checkLinearConsistency(obj,Mesh)
            % Test for Linear Consistency of shape functions
            test=0;
            hw=waitbar(0/Mesh.noElements,'Checking Linear Consistency at Int Points');
            for j=1:Mesh.noElements
                waitbar(j/Mesh.noElements,hw,'Checking Linear Consistency at Int Points');
                for g=1:Mesh.orderInt
                    xTest=Mesh.Elements(j).getIntCord(g); xTest=[xTest(1);xTest(2)]; % Point to Test Consistency
                    ExactResult=xTest(1)+xTest(2);
                    approximateSolution=0;
                    for a=1:obj.numberOfNodes
                        if xTest(1)-obj.Nodes(a).cordinates(1)<=obj.Nodes(a).a && xTest(2)-obj.Nodes(a).cordinates(2)<=obj.Nodes(a).a
                            approximateSolution=approximateSolution+sum(obj.Nodes(a).sF.getValue(xTest)*obj.Nodes(a).cordinates);
                        end
                    end
                    if abs(approximateSolution-ExactResult)>test
                        test=abs(approximateSolution-ExactResult);
                    end
                end
            end
            close(hw)
        end 
        
        function Area=checkArea(obj)
           % Method Sums the area of Voroni Cells
           Area=0;
           for i=1:obj.numberOfNodes
               Area=obj.Nodes(i).Voron.area+Area;
           end
        end
    end
end


