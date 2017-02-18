classdef Voroni < handle
    %VORONI Class containts information on a Voroni Cell (to be associated
    %with an RKNODE Instance)
    
    properties
        VoroniCords; % Counter clockwise description of the Voroni Cell around the node
        normalPoint; % Normal Calculation at midpoints of cells *outward* (the point, used for plotting)
        normal; % Vector describing the normals of the voroni cell
        centre; % Center Point of Voroni Cell
        midPoint; % Midpoint Calculation for cell boundary
        numberOfNodes; % Number of Nodes in the Voroni Cell
        area; % Area of the Cell
        perm; % perimeter of the voroni Cell (total)
        permRatio; % Ratio of perimeter to total perimeter for each side
        radius; % Returns the Maximum radius of the Voroni Cell
        lengths; % Length of Each side
    end
    
    methods 
        % Constructor
        function obj=Voroni(VoroniCords,centre)
           obj.VoroniCords=VoroniCords;
           obj.centre=centre;
           obj.numberOfNodes=size(VoroniCords,1);
           obj.midPoint=zeros(obj.numberOfNodes,2);
           obj.normal=zeros(obj.numberOfNodes,2);
           obj.midPoint=zeros(obj.numberOfNodes,2);
           obj.setNormals()
           obj.setArea();
           obj.setPerimeter();
           obj.setPermRatio();
           obj.setRadius();
        end
        % Determine Normal at MidPoints:
        function setNormals(obj)
            % Method sets the normals of each edge of the Voroni Cell
            Points=obj.VoroniCords;
            Points=[Points;Points(1,:)];   
            for j=1:size(Points,1)-1
                dx=Points(j+1,1)-Points(j,1);
                dy=Points(j+1,2)-Points(j,2);
                obj.midPoint(j,:)=[dx*0.5+Points(j,1),dy*0.5+Points(j,2)];
                r=sqrt(dx^2+dy^2);
                obj.normalPoint(j,:)=obj.midPoint(j,:)+[dy/r, -dx/r];
                obj.normal(j,:)=[dy/r, -dx/r];
            end
        end  
        function setArea(obj)
            % Method sets the area of Voroni Cell
            XY=obj.VoroniCords; XY=[XY;obj.VoroniCords(1,:)];
            area=0;
            for i=1:obj.numberOfNodes
               area=det([XY(i,1) XY(i+1,1);XY(i,2) XY(i+1,2)])+area;
            end
            obj.area=0.5*area;
        end
        
        function setPerimeter(obj)
            % Method determines the perimeter of each Voroni Cell
            Points=obj.VoroniCords;
            Points=[Points;Points(1,:)];   
            obj.perm=0;
            for j=1:size(Points,1)-1
                obj.perm=obj.perm+norm(Points(j,:)-Points(j+1,:));
            end
        end
        
        function setPermRatio(obj)
            % Method determines the fraction of perimeter that each edge
            % belongs to
           obj.permRatio=zeros(obj.numberOfNodes,1);
           obj.lengths=zeros(obj.numberOfNodes,1);
           Points=obj.VoroniCords;
           Points=[Points;Points(1,:)]; 
           for i=1:length(obj.permRatio)
              obj.lengths(i)=norm(Points(i,:)-Points(i+1,:));
              obj.permRatio(i)=obj.lengths(i)/obj.perm;
           end
        end
        
        function setRadius(obj)
           % Method sets the radius (largest edge Midpoint to center distance)
           obj.radius=0;
           for i=1:size(obj.VoroniCords,1)
              rTest=norm(obj.midPoint(i,:)'-obj.centre);
              if rTest>obj.radius
                 obj.radius=rTest; 
              end
           end
        end
        
    end
    
end

