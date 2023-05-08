clear
close all
clc

% Description: 2D finite element method to solve the Magnetostatic
% equations. The output is the magnetic field in space, for given border
% conditions and materials present.

% How to use: Load the mesh desired, establish the parameters 'mu', the 
% border conditions 'g' and 'gam' and run.

% /!\ Important: Currently, density current is not supported

%% Load mesh
load('retangMesh.mat')

t = t(1:3,:);       % mesh (set of triangles)
bnod = find(p(3,:)~=0); % border nodes
nb = max(labels);   % number of borders
nt = length(t);     % number of triangles

%% Define mu in space
mu0 = pi*4e-7;
mu = zeros(nv,1);
mur = 3;    % relative permeability

% Get the object borders
obj_bord = find(p(3,:)==5 |p(3,:)== 6 |p(3,:)==7 |p(3,:)==8);

px = p(1,obj_bord);
py = p(2,obj_bord);

% order the nodes [px,py] by angle to the centroid
% This is for matlab to construct the polygon
meanx = mean(px);
meany = mean(py);
angles = atan2( (py-meany),(px-meanx)); clear meanx meany

[~, sortIndices] = sort(angles);
px = px(sortIndices);
py = py(sortIndices); clear sortIndices angles

pgon = polyshape(px,py);

for nd = 1:nv
    in = inpolygon(p(1,nd),p(2,nd),pgon.Vertices(:,1),pgon.Vertices(:,2));
    if in
        mu(nd) = mu0*mur;
    else
        mu(nd) = mu0;
    end
end

clc
%% Borders and conditions

% Density current
jcurr = zeros(nv,1); 
% jcurrMax = 0;
% for nd = 1:nv
%     [in,on] = inpolygon(p(1,nd),p(2,nd),pgon.Vertices(:,1),pgon.Vertices(:,2));
%     if in || on
%         jcurr(nd) = jcurrMax;
%     end
% end

% Border conditions
g = zeros(1,nb);
gam = g;
g(1) = -1; g(3) = 1;
gam(1) = 0.001; gam(3) = 0.001;

%% (M+R)c = b+j

% j
j = zeros(nv,1);
% for nd = 1:nv
%     tr = t(:,find_triang(nd,t));
% 
%     area = 0;
%     for it = 1:length(tr(1,:))
%         area = area + t_area(tr(:,it),p);
%     end
%     j(nd) = area*jcurr(nd)/3;
% end
     
% M
M = zeros(nv);
for nd1 = 1:nv
    k1 = find_triang(nd1,t);
    nds = unique(t(:,k1));
    for nd2 = 1:nv
        if ismember(nd2,nds)
            k2 = find_triang(nd2,t);
            k = intersect(k1,k2);
            itg = 0;
            for ik = 1:length(k)
                [~,b1,c1,area] = abc(p,t,nd1,k(ik));
                [~,b2,c2,~] = abc(p,t,nd2,k(ik));
                itg = itg + (b1*b2 + c1*c2)*area/(sum(mu(t(:,k(ik))))/3);
            end
            M(nd1,nd2) = itg;
        else
            M(nd1,nd2) = 0;
        end

    end
end

% b
b = zeros(nv,1);

for nd = bnod
    k = find_triang(nd,t);
    nds = unique(t(:,k));
    nds = nds(p(3,nds) == p(3,nd)); % only those on the same border of nd
    nds(nds==nd) = [];

    clear l
    for ind = 1:length(nds)
        l(ind) = sqrt( (p(1,nds(ind))-p(1,nd)).^2 + (p(2,nds(ind))-p(2,nd)).^2 );
    end

    b(nd) = 0.5*sum(l.*g(p(3,nds))./mu(nds)');
end

% R
R = zeros(nv);

for nd1 = bnod

    k1 = find_triang(nd1,t); % all elements of phi i

    nds = unique(t(:,k1));    % nodes of the elements
    nds = nds(p(3,nds) == p(3,nd1)); % only nds of the same border of nd1

    for nd2 = bnod
        if ismember(nd2,nds)
            R(nd1,nd2) = abs(bordIntegral(nd1,nd2,nds,p))*(sum(gam(p(3,nds)))/numel(nds))/(sum(mu(p(3,nds)))/numel(nds));
        end
    end
end

%% Linesolve
coef = linsolve(M+R,j+b);

%% Vector Potential
% Centroid of each element
xs = zeros(nt,1);
ys = xs;

% Pre alocate memory for the magnetic field
Bx = zeros(nt,1);
By = Bx;

% For each element, calculate Bx and By
for k = 1:nt
    nds = t(:,k);

    a = zeros(length(nds),1); b = a; c = a;
    for ind = 1:length(nds)
        nd = nds(ind);
        [at,bt,ct,~] = abc(p,t,nd,k);
        a(ind) = at;
        b(ind) = bt;
        c(ind) = ct;
    end

    a = sum(a.*coef(nds));
    b = sum(b.*coef(nds));
    c = sum(c.*coef(nds));

    Bx(k) = c;
    By(k) = -b;

    % Centroid of the element
    xs(k) = mean(p(1,nds));
    ys(k) = mean(p(2,nds));
end

%% Plot Magnetic Field
figure
hold on
plt = triplot(t',p(1,:),p(2,:),'g');
plt.Color = [plt.Color 0.2];

q = quiver(xs,ys,Bx,By,'b');
q.MaxHeadSize = 0.005;

plot(pgon,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

xlabel('x'); ylabel('y')
% title("Empty Halbach - "+num2str(max(sqrt(Bx.^2 + By.^2)))+" T")
axis('square')

%% The Magnetic material
Inside_elements = zeros(nt,1);
n_elements = 0; % number of elements that make the material
% Store the element index
for k = 1:nt
    in = inpolygon(xs(k),ys(k),pgon.Vertices(:,1),pgon.Vertices(:,2));
    if in
        n_elements = n_elements + 1;
        Inside_elements(n_elements) = k;
    end
end
Inside_elements(Inside_elements==0) = [];

if not(isempty(Inside_elements))
    disp("|B| = "+num2str(max(sqrt(Bx(Inside_elements).^2 + By(Inside_elements).^2)))+" T")
else
    disp("|B| = "+num2str(max(sqrt(Bx.^2 + By.^2)))+" T")
end

%% Functions
function plotri(p,t,i,style)
if nargin < 4
    style = 'b.-';
end
plot([p(1,t(1,i)),p(1,t(2,i)),p(1,t(3,i)),p(1,t(1,i))],[p(2,t(1,i)),p(2,t(2,i)),p(2,t(3,i)),p(2,t(1,i))],style)
end

function t_ind = find_triang(nd,t)
% Objective: Given a node, find all triangle with that node

t_ind = [];

for j = 1:length(t)
    for k = 1:3
        if nd == t(k,j)
            t_ind(end+1) = j;
        end
    end
end

end

function area = t_area(t,p)
% Objective: Given a triangle, calculate the area

for i = 1:length(t)
    pt(:,i) = p(:,t(i));
end

x = pt(1,:); y = pt(2,:);
area = polyarea(x,y);
end

function [a,b,c,area] = abc(p,t,nd,k)
x = p(1,t(t(:,k)~=nd,k));
y = p(2,t(t(:,k)~=nd,k));

area = polyarea([x,p(1,nd)],[y,p(2,nd)]);
a = (x(1)*y(2)-x(2)*y(1))/(2*area);
b = (y(1)-y(2))/(2*area);
c = (x(2)-x(1))/(2*area);
if (a + b*p(1,nd) + c*p(2,nd)) < 0
    a = -a;
    b = -b;
    c = -c;
end

end

function [x,y,m,phi1,phi2] = egde(nd1,nd2,nds,p)

% Objective: Given two nodes, create the function phi i and phi j along
% the edge that connect them.

%% If nd1 ~= nd2
if nd1 ~= nd2
    x = linspace(p(1,nd1),p(1,nd2),20);
    y = linspace(p(2,nd1),p(2,nd2),20);

    if x(end)-x(1) ~= 0
        m = (y(end)-y(1))/(x(end)-x(1));
    else
        m = nan;
    end

else
    nds(nds==nd2) = [];

    % Special case: border node is a corner
    if numel(nds) == 1
        x = [linspace(p(1,nds),p(1,nd1),20);zeros(1,20)];
        y = [linspace(p(2,nds),p(2,nd1),20);zeros(1,20)];

        if x(1,end)-x(1,1) ~= 0
            m = (y(1,end)-y(1,1))/(x(1,end)-x(1,1));
        else
            m = nan;
        end
        m(1:2) = [m,0];
        % Simpler Case: 2 neighbours
    else
        x = [linspace(p(1,nds(1)),p(1,nd1),20);linspace(p(1,nd1),p(1,nds(2)),20)];
        y = [linspace(p(2,nds(1)),p(2,nd1),20);linspace(p(2,nd1),p(2,nds(2)),20)];

        % m1
        if x(1,end)-x(1,1) ~= 0
            m(1) = (y(1,end)-y(1,1))/(x(1,end)-x(1,1));
        else
            m(1) = nan;
        end

        % m1
        if x(2,end)-x(2,1) ~= 0
            m(2) = (y(2,end)-y(2,1))/(x(2,end)-x(2,1));
        else
            m(2) = nan;
        end

    end
end

phi1 = linspace(0,1,20);
phi2 = linspace(1,0,20);

end

function itg = bordIntegral(nd1,nd2,nds,p)

if nd1==nd2
    [x,y,m,phi1,phi2] = egde(nd1,nd2,nds,p);

    % m1
    if isnan(m(1))
        itg = trapz(y(1,:),phi1.*phi1);
    else
        itg = trapz(x(1,:),phi1.*phi1.*sqrt(1+m(1)^2));
    end

    % m2
    if isnan(m(2))
        itg = itg + trapz(y(2,:),phi2.*phi2);
    else
        itg = itg + trapz(x(2,:),phi2.*phi2.*sqrt(1+m(2)^2));
    end

elseif ismember(nd2,nds)
    [x,y,m,phi1,phi2] = egde(nd1,nd2,nds,p);

    if isnan(m)
        itg = trapz(y,phi1.*phi2);
    else
        itg = trapz(x,phi1.*phi2.*sqrt(1+m^2));
    end

end

end
