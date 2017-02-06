%% Create Rectangular Mesh for Quadrature
%Set up Mesh Geometry:
[NN,NEL,X,Y] = GridRectangle(1,1,10,10);
Mesh=QuadMesh(NN,NEL,1);

% Set up Natural Boundary:
b2=[1-eps,1+eps,-eps,1+eps,[0,1]];
BN=Boundary(NN,b2); BN(BN==-Inf)=0;

%% Plot Mesh:
MeshPlot.plotOriginal(Mesh);

%% Generate Uniform Point Cloud:
[NNp] = GridRectangle(1,1,4,4);
%Set up Essential Boundary:
b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b2=[-eps,eps,0-eps,0+eps,[0,-Inf]];
b3=[-eps,eps,1-eps,1+eps,[0,-Inf]];
[Verts,Cells,XY]=VoronoiLimit(NNp(:,2),NNp(:,3),[0,0;1,0;1,1;0,1]); % Create Bounded Voroni Diagram
BE=Boundary([(1:size(XY,1))',XY,reshape((1:size(XY,1)*2),2,[])'],b1,b2,b3);
xNodes=XY(:,1);
yNodes=XY(:,2);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,1,Verts,Cells);
PointCloud.plotCloud()
PointCloud.plotVoroni(0.01)

%% Define Random Point Cloud Object:
Points=Mesh.createPoints(50,[0,0;0,0.5;0,1;1,1;1,0]); % Create Initial Points
[Verts,Cells,XY]=VoronoiLimit(Points(:,1),Points(:,2),[0,0;1,0;1,1;0,1]); % Create Bounded Voroni Diagram
%Set up Essential Boundary:
b1=[-eps,eps,0.5-eps,0.5+eps,[0,0]];
b2=[-eps,eps,0-eps,0+eps,[0,0]];
b3=[-eps,eps,1-eps,1+eps,[0,0]];
BE=Boundary([(1:size(XY,1))',XY,reshape((1:size(XY,1)*2),2,[])'],b1,b2,b3);
% Pass in Points from the XY data of the Voroni Cells
xNodes=XY(:,1);
yNodes=XY(:,2);
Nodes=[(1:length(xNodes))',xNodes,yNodes,reshape(BE,2,[])'];
PointCloud=Cloud(Nodes,1,2,Verts,Cells);
PointCloud.plotCloud()
PointCloud.plotVoroni(0.01)

%% Perform Checks
PointCloud.checkArea()
PointCloud.checkPU()
%PointCloud.checkNU(Mesh)
%PointCloud.checkLinearConsistency(Mesh)

%% Build Arrays:
C=Constit(1E4,0,'Plane Stress').C;
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
PointCloud.returnInterpolatedU([1;0])

%%
PointCloud.plotVoroniDeformed(1)
