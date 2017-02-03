classdef Voroni < handle
    %VORONI Class containts information on a Voroni Cell (to be associated
    %with an RKNODE Instance)
    
    properties
        VoroniCords; % Counter clockwise description of the Voroni Cell around the node
        normalPoint; % Normal Calculation at midpoints of cells *outward* (the point, used for plotting)
        normal; % Vector describing the normals of the voroni cell
        midPoint; % Midpoint Calculation for cell boundary
        numberOfNodes; % Number of Nodes in the Voroni Cell
        area; % Area of the Cell
    end
    
    methods 
        % Constructor
        function obj=Voroni(VoroniCords)
           obj.VoroniCords=VoroniCords;
           obj.numberOfNodes=size(VoroniCords,1);
           obj.midPoint=zeros(obj.numberOfNodes,2);
           obj.normal=zeros(obj.numberOfNodes,2);
           obj.midPoint=zeros(obj.numberOfNodes,2);
           obj.setNormals()
           obj.setArea();
        end
        % Determine Normal at MidPoints:
        function setNormals(obj)
            Points=obj.VoroniCords;
            Points=[Points;Points(1,:)];   
            for j=1:size(Points,1)-1
                dx=Points(j+1,1)-Points(j,1);
                dy=Points(j+1,2)-Points(j,2);
                obj.midPoint(j,:)=[dx/2+Points(j,1),dy/2+Points(j,2)];
                r=sqrt(dx^2+dy^2);
                obj.normalPoint(j,:)=obj.midPoint(j,:)+[dy/r, -dx/r];
                obj.normal(j,:)=[dy/r, -dx/r];
            end
        end  
        function setArea(obj)
            XY=obj.VoroniCords; XY=[XY;obj.VoroniCords(1,:)];
            area=0;
            for i=1:obj.numberOfNodes
               area=det([XY(i,1) XY(i+1,1);XY(i,2) XY(i+1,2)])+area;
            end
            obj.area=0.5*area;
        end
    end
    
end

