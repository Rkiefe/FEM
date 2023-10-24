clear
close all
clc

% Basic physical properties
mu0 = pi*4e-7;  % permeability of vaccuum

mur = 1.05;     % relative permeability

%% Mesh
% Mesh file name
mesh.name = "mesh";

% Load mesh
load("../Mesh/"+mesh.name+".mat")

% change border name
borders = b; clear b ans

% number of borders
nb = max(labels);

% border nodes
bnod = find(p(3,:)~=0);

% change the triangle matrix to only have the nodes
t = t(1:3,:);

% xs and ys || centroid of each element k
xs = zeros(nt,1);
ys = xs;

for k = 1:nt
    nds = t(:,k);

    % Centroid of the element
    xs(k) = mean(p(1,nds));
    ys(k) = mean(p(2,nds));
end

%% Object - Magnetizable object
obj_bord = find(p(3,:)== 5 |p(3,:)== 6 |p(3,:)== 7 |p(3,:)== 8);

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

% Elements that belong to the object:
Inside_elements = zeros(nt,1);
n_elements = 0; % number of elements that make the material (will increase)
% Store the element index
for k = 1:nt
    in = inpolygon(xs(k),ys(k),pgon.Vertices(:,1),pgon.Vertices(:,2));
    if in
        n_elements = n_elements + 1;
        Inside_elements(n_elements) = k;
    end
end
Inside_elements(n_elements+1:end) = []; clear n_elements


%% Simulation Parameters

% Border conditions
g = zeros(1,nb);

Hext = 1.1/mu0; % strength of the magnetic field in A/m --> Tesla / mu0

g(6) = -Hext;
g(8) = Hext;


%% Start FEM-Solver
mu = zeros(nt,1) + mu0;

% Get chi and mu on each element
mu(Inside_elements) = mu0*mur;

% A matrix
A = zeros(nv);
for nd1 = 1:nv
    k1 = find_triang(nd1,t);

    nds = unique(t(:,k1));
    for nodeInd = 1:numel(nds)
        nd2 = nds(nodeInd);

        k2 = find_triang(nd2,t);
        k = intersect(k1,k2);

        for ik = 1:numel(k)
            [~,bi,ci,area] = abc(p,t,nd1,k(ik));
            [~,bj,cj,~] = abc(p,t,nd2,k(ik));

            A(nd1,nd2) = A(nd1,nd2) + (bi*bj+ci*cj)*area*mu(k(ik));
        end
    end
end

% q vector
q = zeros(nv,1);
for ind = 1:numel(bnod)
    nd = bnod(ind);

    % find adjacent nodes to 'nd'
    lst = unique([find(borders(1,:)==nd),find(borders(2,:)==nd)]);
    nds = unique(borders(1:2,lst)); nds(nds==nd) = [];

    % border nodes and indices
    g_ind = [nds,borders(3,lst)'];

    % calculate the integral
    q(nd) = qVec(nd,p,g*mu0,g_ind);
end

% Lagrange multiplier

% 1st - Calculate C
C = zeros(1,nv);
for nd = 1:nv

    k_list = find_triang(nd,t); % list of elements with node nd
    aux = 0;
    for ki = 1:numel(k_list)
        k = k_list(ki);

        aux = aux + polyarea(p(1,t(:,k)) , p(2,t(:,k)))/3;
    end
    C(nd) = aux;
end

mat = zeros(length(A)+1);
mat(1:length(A),1:length(A)) = A;
mat(end,1:length(C)) = C;
mat(1:length(C),end) = C';

r = [-q;0];
x = mat\r;      % lagrange multiplier

coef = x(1:end-1);  % FEM result

% New field H
H = zeros(nt,2);
for k = 1:nt
    nds = t(:,k); % all nodes of that element

    % Sum the contributions
    for ind = 1:numel(nds)
        nd = nds(ind);
        % obtain the element parameters
        [~,bi,ci,~] = abc(p,t,nd,k);

        H(k,1) = H(k,1) - coef(nd)*bi;
        H(k,2) = H(k,2) - coef(nd)*ci;
    end
end

% Chi
chi = mu./mu0 -1;

% M
M = chi.*H;

% B
B = mu0.*(H+M);

%% Plot
fig = figure;
tiledlayout(1,3);

plotResult(t,p,pgon,xs,ys,B(:,1),B(:,2),"T","B",fig)
plotResult(t,p,pgon,xs,ys,H(:,1),H(:,2),"A/m","H",fig)
plotResult(t,p,pgon,xs,ys,M(:,1),M(:,2),"A/m","M",fig)

%% Compare with FEMM
% close all

% \\\ New
% Select data along line of interest
ind = [];
for k = 1:nt
    if abs(ys(k) - 0) < 0.25 && ( xs(k) < 3 && xs(k) > -3 )
        ind = [ind,k];
    end
end

x = xs(ind);
yB = sqrt(sum(B(ind,:).^2,2)); % Sort B
yH = sqrt(sum(H(ind,:).^2,2)); % Sort H

[x,i] = sort(x);
yB = yB(i);
yH = yH(i);

F_B = readmatrix('../FEMM_B');
F_H = readmatrix('../FEMM_H');

u = linspace(-3,3,1e3);
aux = interp1(x,yB,u); yB = aux;
aux = interp1(x,yH,u); yH = aux;

aux = interp1(F_B(:,1)-3,F_B(:,2),u); F_B = aux;
aux = interp1(F_H(:,1)-3,F_H(:,2),u); F_H = aux;

figure
subplot(1,2,1)
plot(u,yB,'.',u,F_B,'.')

title('B')

subplot(1,2,2)
plot(u,yH,'.',u,F_H,'.')

title('H')

%% Functions
function [chi,H_refined] = getChi(T_sample)
prof = load('../Gd_simulado/workspace_Gd_logH_2T');
mu0 = pi*4e-7;
rho = 7.9; % g/cm3

H_T = prof.H_T;
T_K = prof.T_K;
M_minG_smooth_emug = prof.M_minG_smooth_emug;

H_refined = linspace(min(H_T),max(H_T),1e3)./mu0;

indT = find(abs(T_K-T_sample) == min(abs(T_K-T_sample)),1);

M_Am = M_minG_smooth_emug.*rho.*1e3;  % A/m
H_Am = H_T./mu0;

chi = M_Am(:,indT)./H_Am;
chi = interp1(H_Am,chi,H_refined);
end

function itg = qVec(nd,p,g,g_ind)
nds = g_ind(:,1);
b_idx = g_ind(:,2);

x = [linspace(p(1,nds(1)),p(1,nd),20);linspace(p(1,nd),p(1,nds(2)),20)];
y = [linspace(p(2,nds(1)),p(2,nd),20);linspace(p(2,nd),p(2,nds(2)),20)];

% m1
if x(1,end)-x(1,1) ~= 0
    m(1) = (y(1,end)-y(1,1))/(x(1,end)-x(1,1));
else
    m(1) = nan;
end

% m2
if x(2,end)-x(2,1) ~= 0
    m(2) = (y(2,end)-y(2,1))/(x(2,end)-x(2,1));
else
    m(2) = nan;
end
phi1 = linspace(0,1,20);
phi2 = linspace(1,0,20);

% Integral:

% m1
if isnan(m(1))
    itg = g(b_idx(1))*abs(trapz(y(1,:),phi1));
else
    itg = g(b_idx(1))*abs(trapz(x(1,:),phi1)*sqrt(1+m(1)^2));
end

% m2
if isnan(m(2))
    itg = itg + g(b_idx(2))*abs(trapz(y(2,:),phi2));
else
    itg = itg + g(b_idx(2))*abs(trapz(x(2,:),phi2)*sqrt(1+m(2)^2));
end

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

function plotResult(t,p,pgon,xs,ys,Bx,By,units,field,fig)
figure(fig)
tl = nexttile;

plt = triplot(t',p(1,:),p(2,:),'g'); hold on
plt.Color = [plt.Color 0.2];
q = quiver(xs,ys,Bx,By,'b');
q.MaxHeadSize = 0.005;
plot(pgon,'FaceAlpha',0,'EdgeColor','r','LineWidth',1)

Bmax = max(sqrt(Bx.^2 + By.^2));

title(field+"_{max} = "+num2str(Bmax)+" "+units)
xlabel('x'); ylabel('y')
end