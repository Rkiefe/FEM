function out = FEM3D_linear(model,p,t,surfaceT,surface2element,Hext,direction)

if nargin < 7
    direction = "vertical";
end

msh = model.Mesh;

% Volume of each element of the mesh
[~,VE] = volume(msh); % Volume of geometry and volume of each element

nv = length(p); % number of nodes
nt = length(t); % number of tetrahedra

surfaceNodeInd = find(p(4,:)~=0);

%% Center coordinates of each element

xc = zeros(nt,1);
yc = zeros(nt,1);
zc = zeros(nt,1);

for elementIndex = 1:nt

    % Nodes of that element
    nodesElement = t(:,elementIndex);

    xc(elementIndex) = mean(p(1,nodesElement));
    yc(elementIndex) = mean(p(2,nodesElement));
    zc(elementIndex) = mean(p(3,nodesElement));
end

pc = [xc,yc,zc];

out.centerCoord = pc;

% Elements that belong inside the object:
InsideElements = findElements(msh,"region","Cell",2); % Elements from Cell 2

out.InsideElements = InsideElements;
%% Border Conditions and system parameters

mu0 = pi*4e-7;  % Vaccum permeability
mur = 10;       % Relative permeability

g = zeros(model.Geometry.NumFaces,1);
gam  = zeros(model.Geometry.NumFaces,1) + 1e-9; % + 1e-9

if direction == "horizontal"
    g(5) = -Hext;
    gam(3) = 1e6;
else
    g(1) = -Hext;
    gam(2) = 1e6;
end

% Permeability in space
mu = zeros(nt,1) + mu0;
mu(InsideElements) = mu0*mur;

%% Finite Element matrices

% Aij
A = zeros(nv);
for nd1 = 1:nv

    % All elements with node nd1
    nd1List = find_elements(nd1,t);

    % Only the neighbouring nodes give Aij ~= 0
    possibleNodes = t(:,nd1List); possibleNodes = unique(possibleNodes(:));

    for ind = 1:numel(possibleNodes)
        nd2 = possibleNodes(ind);

        % All elements with node nd2
        nd2List = find_elements(nd2,t);

        % All shared elements between nd1 and nd2
        sharedElements = intersect(nd1List,nd2List);

        % For each shared element, calculate the integral
        itg = 0;
        for elInd = 1:numel(sharedElements)

            % Current element index
            currElement = sharedElements(elInd);

            % a b c d parameters
            [~,bi,ci,di] = abcd(p,t(:,currElement),nd1);
            [~,bj,cj,dj] = abcd(p,t(:,currElement),nd2);

            % Volume of the tetrahedral
            Vk = VE(currElement);

            itg = itg + mu(currElement)*(bi*bj + ci*cj + di*dj)*Vk;
        end
        A(nd1,nd2) = itg;
    end
end

% qj
q = zeros(nv,1);

% Only need to scan the nodes belonging to a surface
for ind = 1:numel(surfaceNodeInd)

    % Current Node
    node = surfaceNodeInd(ind);

    % Find all surface triangles with current node
    triangleIndex = find_surfaceTriangles(node,surfaceT);

    % Integrate for each surface triangle
    itg = 0;
    for it = 1:numel(triangleIndex)
        % Current Triangle
        currTriangle = triangleIndex(it);

        % Element associated with that surface triangle
        currElement = surface2element(currTriangle);

        % Obtain psi parameters a b c d
        [ai,bi,ci,di] = abcd(p,t(:,currElement),node);

        % Triangle vertices coordinates
        xt = p(1,surfaceT(1:3,currTriangle));
        yt = p(2,surfaceT(1:3,currTriangle));
        zt = p(3,surfaceT(1:3,currTriangle));

        % Centroid of the triangle
        xcTri = mean(xt);
        ycTri = mean(yt);
        zcTri = mean(zt);

        % Area of the current surface triangle
        area = areaTriangle(xt,yt,zt);

        % Border of the current triangle
        bCurrent = surfaceT(4,currTriangle);

        % Calculate psi on the centroid of the triangle
        psi = ai + bi*xcTri + ci*ycTri + di*zcTri;

        itg = itg + mu0*psi*g(bCurrent)*area;
    end
    q(node) = itg;
end

% Rij
R = zeros(nv);
% Only need to scan the nodes belonging to a surface
for ind = 1:numel(surfaceNodeInd)

    % Current nd1
    nd1 = surfaceNodeInd(ind);

    % Find all surface triangles with node nd1
    trianglesNd1 = find_surfaceTriangles(nd1,surfaceT);

    % All nodes of the surface triangles that contain nd1
    nodes = surfaceT(1:3,trianglesNd1); nodes = unique(nodes(:));

    for jnd = 1:numel(nodes)

        % Current node nd2
        nd2 = nodes(jnd);

        % Surface triangles that have nd2
        trianglesNd2 = find_surfaceTriangles(nd2,surfaceT);

        % Shared triangles
        triangleShared = intersect(trianglesNd1,trianglesNd2);

        % For each shared triangle, calculate the integral
        itg = 0;
        for it = 1:numel(triangleShared)
            % Current Triangle
            currTriangle = triangleShared(it);

            % Corresponding element
            currElement = surface2element(currTriangle);

            % Obtain psi parameters a b c d
            [ai,bi,ci,di] = abcd(p,t(:,currElement),nd1);
            [aj,bj,cj,dj] = abcd(p,t(:,currElement),nd2);

            % Triangle vertices coordinates
            xt = p(1,surfaceT(1:3,currTriangle));
            yt = p(2,surfaceT(1:3,currTriangle));
            zt = p(3,surfaceT(1:3,currTriangle));

            % Centroid of the triangle
            xcTri = mean(xt);
            ycTri = mean(yt);
            zcTri = mean(zt);

            % Area of the current surface triangle
            area = areaTriangle(xt,yt,zt);

            % Border of the current triangle
            bCurrent = surfaceT(4,currTriangle);

            % Calculate psi on the centroid of the triangle
            psiI = ai + bi*xcTri + ci*ycTri + di*zcTri;
            psiJ = aj + bj*xcTri + cj*ycTri + dj*zcTri;
            
            itg = itg + mu0*gam(bCurrent)*psiI*psiJ*area;

        end
        R(nd1,nd2) = itg;
    end


end


%% Resulting Magnetic field
coef = (A-R)\(-q);

% New field H
H = zeros(nt,3);
for elInd = 1:nt

    elNodes = t(:,elInd); % all nodes of that element

    % Sum the contributions
    for ind = 1:numel(elNodes)
        nd = elNodes(ind);
        
        % obtain the element parameters
        [~,bi,ci,di] = abcd(p,elNodes,nd);

        H(elInd,1) = H(elInd,1) - coef(nd)*bi;
        H(elInd,2) = H(elInd,2) - coef(nd)*ci;
        H(elInd,3) = H(elInd,3) - coef(nd)*di;
    end
end

% B field
B = zeros(nt,3);
B(:,1) = H(:,1).*mu;
B(:,2) = H(:,2).*mu;
B(:,3) = H(:,3).*mu;

out.H = H;
out.B = B;

%% Functions
function A = areaTriangle(xt,yt,zt)

ons = [1 1 1];

A = 0.5*sqrt(det([xt;yt;ons])^2 + det([yt;zt;ons])^2 + det([zt;xt;ons])^2);
end

function tind = find_surfaceTriangles(nd,surfaceT)
% Objective: Given a node, find all tetrahedral with that node

tind = [];

for j = 1:length(surfaceT)
    for ndInd = 1:3
        if nd == surfaceT(ndInd,j)
            tind(end+1) = j;
        end
    end
end

end

function tind = find_elements(nd,t)
% Objective: Given a node, find all tetrahedral with that node

tind = [];

for j = 1:length(t)
    for ndInd = 1:4
        if nd == t(ndInd,j)
            tind(end+1) = j;
        end
    end
end

end

function [a,b,c,d] = abcd(p,nodes,nd)

nodes(nodes==nd) = [];

x = [p(1,nd),p(1,nodes)]';
y = [p(2,nd),p(2,nodes)]';
z = [p(3,nd),p(3,nodes)]';

psi = [1;0;0;0];

M = [ones(4,1),x,y,z];

aux = M\psi;

a = aux(1);
b = aux(2);
c = aux(3);
d = aux(4);

end

end