classdef RKNode < handle
    %NODE is an object for an RKPM Mesh
    
    properties
        sF; % RKPM Shape Function
        weight; % Weight Function
        a; % Dilation Parameter
        order; % order of approximation
        nodeNumber; % Global Node Number
        cordinates; % List of Cordinates for the Node
        BE; % Flag for the essential Boundary Status
        singularKernal=false; % Flag for whether or not Node will use Singular Kernal
        Cloud; % Cloud Data Handle
        u; % Solutions *** Lack of Kronecker Means this isn't the PDE Solution Yet
        Voron; % Object containing the voroni Diagram
        CellDeriv; % Data Structure containing the Cata
    end
    
    methods
        function obj = RKNode(nodeNumber,cordinates,dilation,order,BE,VoroniCords,Cloud)
            % Set Attributes:
            obj.nodeNumber = nodeNumber;
            obj.cordinates=cordinates;
            obj.a=dilation;
            obj.BE=BE;
            obj.order=order;
            obj.weight=[];
            obj.Cloud=Cloud;
            if sum(BE~=-Inf)
                obj.singularKernal=true;
            end
            % Define Soldier Classes: 
            obj.setWeightFunction()
            obj.setShapeFunction()
            obj.setVoroni(VoroniCords)
            
            % Define Derivative Data Structure
            obj.CellDeriv=CellDeriv(obj.Cloud,obj.Voron);
            
            
        end
        % Sets Initial Shape Function
        function setShapeFunction(obj)
            obj.sF=RKShape(obj.cordinates,obj.a,obj.order,obj.Cloud,obj.weight);
        end
        
        % Sets Initial Weight Function
        function setWeightFunction(obj)
            obj.weight=Weight(obj.cordinates,obj.a,obj.singularKernal);
        end
        
        % Set Value of U
        function setU(obj,u)
            obj.u=u;
        end
        
        % Set Value of Voroni
        function setVoroni(obj,VoroniCords)
            obj.Voron=Voroni(VoroniCords,obj.cordinates);
        end
    end
    
end

