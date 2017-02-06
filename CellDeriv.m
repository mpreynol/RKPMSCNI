classdef CellDeriv < handle
    %CELLDERIV stores the information of all other nodes derivatives at the
    %node this instance lives in
    
    properties
        B; % Stored Derivatives of all nodes (nx2)
        Cloud; % Data Handle for the Point Cloud
        Voron; % Data Handle for this Nodes Voroni Cell
        n; % Number of Nodes
        F; % Force Vector
    end
    
    methods
        function obj=CellDeriv(Cloud,Voron)
            obj.Cloud=Cloud;
            obj.Voron=Voron;
            obj.n=obj.Cloud.numberOfNodes;
            obj.B=zeros(obj.n,2);
            obj.F=zeros(obj.n*2,1);
        end
        
        function setData(obj,Q)
                A=obj.Voron.area;
                listPerm=obj.Voron.permRatio;
                cellCentre=obj.Voron.centre;
                cellRadius=obj.Voron.radius;
                n1=[1;0]; n2=[0;1];
                for a=1:obj.n % Loop over Shape Functions (Rows)
                    cordsA=obj.Cloud.Nodes(a).cordinates;
                    if norm(cellCentre-cordsA)<=obj.Cloud.Nodes(a).a+cellRadius % Then Point a may be included
                        b_a1=0; % Initialize a's Component of the B Matrix
                        b_a2=0; % Initialize a's Component of the B Matrix
                        for l=1:obj.Voron.numberOfNodes % Loop through edges
                            midPoint=(obj.Voron.midPoint(l,:))';
                            normal=(obj.Voron.normal(l,:))';
                            va=obj.Cloud.Nodes(a).sF.getValue(midPoint);
                            b_a1=b_a1+1/A*va*(dot(normal,n1))*obj.Voron.lengths(l);
                            b_a2=b_a2+1/A*va*(dot(normal,n2))*obj.Voron.lengths(l);
                            if sum(Q~=0)>0
                                obj.F(2*a-1,1)=obj.F(2*a-1,1)+va*A*listPerm(l)*Q(1);
                                obj.F(2*a,1)=obj.F(2*a,1)+va*A*listPerm(l)*Q(2);
                            end
                        end
                        obj.B(a,:)=[b_a1,b_a2];
                    end
                end
        end   
    end
end

