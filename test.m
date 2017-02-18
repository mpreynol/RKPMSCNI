%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(10,2,10,20);
Mesh=QuadMesh(NN,NEL,3);

% Set up Natural Boundary:
b2=[10-eps,10+eps,-eps,2+eps,[0,-1]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh);

%% Generate Uniform Point Cloud:
[NNp] = GridRectangle(10,2,25,7);
%Set up Essential Boundary:
b1=[-eps,eps,-eps,2+eps,[0,0]];
[Verts,Cells,XY]=VoronoiLimit(NNp(:,2),NNp(:,3)); % Create Bounded Voroni Diagram
BE=ExactBoundary([(1:size(XY,1))',XY,reshape((1:size(XY,1)*2),2,[])'],b1);
xNodes=XY(:,1);
yNodes=XY(:,2);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,1.4,Verts,Cells);
PointCloud.plotCloud()
PointCloud.plotVoroni(0.03)
PointCloud.checkArea()

%% Define Random Point Cloud Object:
Points=Mesh.createPoints(120,[0,0;0,0.3;0,1;0,1.5;0,2;10,0;10,1;10,2;0,0.75;0,1.7;0,1.2;10,0.75;10,1.7;10,1.2]); % Create Initial Points
[Verts,Cells,XY]=VoronoiLimit(Points(:,1),Points(:,2),[0,0;1,0;1,1;0,1]); % Create Bounded Voroni Diagram
%Set up Essential Boundary:
b1=[-eps,eps,-eps,2+eps,[0,0]];
BE=ExactBoundary([(1:size(XY,1))',XY,reshape((1:size(XY,1)*2),2,[])'],b1);
% Pass in Points from the XY data of the Voroni Cells
xNodes=XY(:,1);
yNodes=XY(:,2);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,2,Verts,Cells);
PointCloud.plotCloud()
PointCloud.plotVoroni(0.01)
PointCloud.checkArea()

%% Perform Checks
PointCloud.checkArea()
PointCloud.checkPU()  
%PointCloud.checkNU(Mesh)
%PointCloud.checkLinearConsistency(Mesh)

%% Build Arrays:
C=Constit(21.1E6,0.3,'PlaneStress').C;
Q=[0;0];
[K,Fb]=PointCloud.integrateDomain(C,Q); % Build Arrays from Voroni
Fh=PointCloud.integrateExactBoundary(Mesh,BN,@parabolicStress);
F=Fb+Fh; sum(F)
%% Solve System:
L=BE==-inf; % Indexes of unknown equations
Kr=K(L,L); Br=BE(~L); fr=F(L); KRHS=K(L,~L); RHS=fr-KRHS*Br;
ur=Kr\RHS;
u=PointCloud.reAssembleUnknowns(ur,BE);
% Populate solution back into PointCloud Collection:
PointCloud.parseSolution(u);

%% Plot The Solution:
PointCloud.plotCloud();
PointCloud.plotCloudU(1);

% %% Test Points
% xs=0:0.1:1;
% for q=1:length(xs)
%     ui=PointCloud.returnInterpolatedU([xs(q);1]);
%     plot(xs(q),ui(2),'k.');
%     hold on
% end
% 

%%
PointCloud.returnInterpolatedU([10;0])

%%
PointCloud.plotVoroniDeformed(1)
%% Section Cut
xSample=(0:0.1:10)';
uSample=zeros(length(xSample),1);
YSpot=1;
% Assemble Data Arrays
for w=1:length(xSample)
    A=PointCloud.returnInterpolatedU([xSample(w);YSpot]);
    uSample(w)=A(2);
end
