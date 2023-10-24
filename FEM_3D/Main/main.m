clear
close all
clc

%% 3D Model and mesh

% Make an empty model
model = createpde;

% Make Box of air
g1 = multicuboid(7,7,8);

% Make Gd geometry                  Bhal            Ours
g2 = multicuboid(0.9,2.5,4); % (0.9,2.5,4)   (1.65,1.65,0.08)
g2 = translate(g2,[0 0 2]); % Move it to the center of the air box

% Make the final geometry
gm = addCell(g1,g2);


% Add the geometry to the model
model.Geometry = gm;

% Number of faces in the geometry
numFaces = model.Geometry.NumFaces;

% Target maximum element size for the mesh
hMax = 1;
% Target minimum element size for the mesh
hMin = 0.1;

% Growth of the element size
hGrow = 1.2;

% Generate a mesh for the new model
msh = generateMesh(model,"GeometricOrder","linear","Hmax",hMax,"Hmin",hMin,"Hgrad",hGrow); % ,"Hmin",hMin

p = msh.Nodes;
t = msh.Elements;

nv = length(p); % number of nodes
nt = length(t); % number of elements

n2 = numel(findElements(msh,"region","Cell",2));
n1 = numel(findElements(msh,"region","Cell",1));

disp("Number of elements: "+num2str(nt))
disp("Number of Inside elements: "+num2str(n2))
pause

% /// Plot Object and mesh
mshPlot = figure;
subplot(1,2,1)
pdemesh(model,FaceAlpha=0.1);

% Plot face labels
subplot(1,2,2)
pdegplot(model,"FaceLabels","on","EdgeLabels","off","CellLabels","on","FaceAlpha",0.1);

%%

%  /// Get the surface nodes and triangles
% /!\ Note: Nodes can belong to more than one face

p = [p;zeros(1,nv)]; % add a row for node border index

% Create a border cell to store all the faces each node belongs to
borderCell = cell(nv,1);

% Find nodes on Face f
for f = 1:model.Geometry.NumFaces

    % Find the nodes that belong to face "f"
    Nf = findNodes(msh,"region","Face",f);

    for ind = 1:numel(Nf)
        if ~isempty(borderCell{Nf(ind)})
            borderCell{Nf(ind)} = [f,borderCell{Nf(ind)}];
        else
            borderCell{Nf(ind)} = f; % Store the node border index
        end
    end

    p(4,Nf) = f;
end

% Nodes that belong to a surface
surfaceNodeInd = find(p(4,:)~=0);

% /// From the nodes, obtain the surface triangles

% Matrix with surface nodes that compose a triangle, and the border index
surfaceT = []; % [nd1,nd2,nd3,border]'

% Map each surface triangle to element
surface2element = [];

% /// For each element, evaluate the border index of each of the nodes: p(4,nodes)
for elInd = 1:nt

    % Nodes
    elNodes = t(:,elInd);

    % Evaluate the border index of each node
    borderArray = borderCell(elNodes);

    % Select only the border index that belong to a surface border
    [ndBord,localIndex] = belong2Surface(borderArray,elNodes);

    % If less than 3 are in a surface border, move on
    if numel(ndBord) < 3
        continue
    end

    % Check the number of nodes that belong to each face
    counts = bordCount(borderArray,numFaces);

    % Check what surface has a triangle
    s = [];
    for ind = 1:numel(counts)
        if counts(ind) > 2 % It means there are 3 nodes in one single surface - a triangle
            s = [s,ind];
        end
    end

    % For each surface, scan what nodes belong to that surface
    for ind = 1:numel(s)
        sCurr = s(ind);

        % New Triangle
        tr = [];

        % Scan each node...
        for jnd = 1:numel(ndBord)
            i = localIndex(jnd); % Node index in the border Array

            % Check if that node belongs to current surface s
            aux = intersect(sCurr,borderArray{i});

            % Store the node to make the triangle
            if ~isempty(aux)
                tr = [tr;ndBord(jnd)];
            end

        end

        % add the surface from which the triangles belongs too
        tr = [tr;sCurr];

        % Save the new triangle
        surfaceT = [surfaceT,tr];

        % Map the new surface triangle to its element
        surface2element(end+1) = elInd;
    end

end

% Remove the repeated surface Triangles
[~,idx]=unique(sort(surfaceT(1:3,:)',2),'rows','stable');
surfaceT = surfaceT(:,idx);

%% FEM

mu0 = pi*4e-7;  % Vaccum permeability
Hext = 1.1/mu0; % external magnetic field

% direction = "horizontal";
% FEM3D(model,p,t,surfaceT,surface2element,Hext,TSample,direction);

direction = "vertical";
out = FEM3D_linear(model,p,t,surfaceT,surface2element,Hext,direction);

%% Plot vector field
H = out.H;
B = out.B;
pc = out.pc;
InsideElements = out.InsideElements;

Hplot = figure;
subplot(1,2,1)
plotVectorField(H,surfaceT,p,pc,"H",[-2 2 -2 2 1 7])

subplot(1,2,2)
plotVectorField(B,surfaceT,p,pc,"B",[-2 2 -2 2 1 7])

%% Plot field Intensity

modH = sqrt(sum(H.^2,2)).*mu0;

figure(mshPlot)
subplot(1,2,2)
hold on
scatter3(pc(InsideElements,1),pc(InsideElements,2),pc(InsideElements,3),30,modH(InsideElements), "filled")
colorbar


%% Functions
function plotVectorField(H,surfaceT,p,pc,text,axi)
if nargin < 5 % No field label
    text = "H"; % Default to H
end

% Plot Vector field
q_plt = quiver3(pc(:,1),pc(:,2),pc(:,3),H(:,1),H(:,2),H(:,3));
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
title("Field "+text)

hold on
ax = gca; ax.FontSize = 20;

for it = 1:length(surfaceT)
    if surfaceT(4,it) == 7 || surfaceT(4,it) == 8 || surfaceT(4,it) == 8 || surfaceT(4,it) == 9 || surfaceT(4,it) == 10 || surfaceT(4,it) == 11 || surfaceT(4,it) == 12
        tr = surfaceT(1:3,it);
        patch = fill3(p(1,tr),p(2,tr),p(3,tr),'r');
        patch(1).FaceAlpha = 0.5;
    end
end

if nargin < 6
    % Zoom in on the area of interest
    axis([-2.5 2.5 -2.5 2.5 4 7]);
else
    axis(axi);
end

end

function counts = bordCount(borderArray,numFaces)
bins = 1:numFaces+1;
occur = [];
for ind = 1:numel(borderArray)
    occur = [occur,borderArray{ind}];
end

% Count how many nodes belong to each surface
counts = histcounts(occur,bins);

end

function [ndBord,localIndex] = belong2Surface(borderArray,elNodes)
ndBord = [];
localIndex = [];
for ind = 1:numel(borderArray)
    if ~isempty(borderArray{ind})
        ndBord = [ndBord,elNodes(ind)];
        localIndex = [localIndex,ind];
    end
end

end
