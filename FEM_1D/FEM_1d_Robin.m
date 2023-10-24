% Author: Rodrigo Kiefe
% Date: 15/02/2023

% Purpose:
% This code simulates the displacement of an elastic bar with Robin
% boundary conditions

% Recommended Reading:
% Mats G.Larson - Theory Implementation and Applications
% Z. Chen - Finite Element Methods and Their applications

clear
close all
clc

%% Setup
% Bar length
L = 1;

% Mesh
nv = 25;    % number of nodes / vertices
nod = linspace(0,L,nv); % noode coordinates

% Force Applied
% f = ones(nv,1);
% f = 1./(1+20.*nod);
f = nod.^3;

kspring = 1e3; % spring constant

g0   = 0;  % displacement of the bar at 0
gL   = 0;  % displacement of the bar at L


%% Boundary and conditions
% Domain \Gamma --- x \in \Gamma
res = 1e-2;
x = 0:res:L;


%% Trial function phi:

phi = trial_fun(x,nod);

dphi = zeros(numel(x),nv);
for in = 1:nv
    dphi(:,in) = gradient(phi(:,in))/res;
end

% (A + R)c = F + r

%% Matrix A
A = zeros(nv);
for i = 1:nv
    for j = 1:nv
        A(i,j) = trapz(x,dphi(:,i).*dphi(:,j));
    end
end

%% Matrix R
R = zeros(nv);
for i = 1:nv
    for j = 1:nv
        R(i,j) = kspring*phi(end,j)*phi(end,i) + kspring*phi(1,j)*phi(1,i);
    end
end

%% F vector

F = zeros(nv,1);
for i = 1:nv
    F(i) = trapz(x,f(i)*phi(:,i));
end

%% r vector

r = zeros(nv,1);
for i = 1:nv
    r(i) = kspring*gL*phi(end,i) + kspring*g0*phi(1,i);
end

%% Solve system of equations:

AR = A+R;
Fr = F+r;
coef = linsolve(AR,Fr);

% solution:
u = coef;

%% Plot
initial = linspace(g0,gL,nv);

subplot(2,1,1)
plot(nod,u,'b.',nod,initial,'go')
xlabel('x'); ylabel('u(x)')
title('1D Elastic Rod')
legend("Displacement", "Starting Position")

subplot(2,1,2)
plot(nod,f,'k-')
xlabel('x'); ylabel('Force Intensity');
title('Force applied')

%% Function
function phi = trial_fun(x,nod)
nv = numel(nod);
phi = zeros(numel(x),nv);

for in = 2:nv-1
    h1 = nod(in)-nod(in-1);
    h2 = nod(in+1)-nod(in);

    k = 1;
    for xc = x
        if nod(in-1) < xc && xc < nod(in)
            phi(k,in) = (xc-nod(in-1))/h1;
        elseif nod(in) < xc && xc < nod(in+1)
            phi(k,in) = (nod(in+1)-xc)/h2;
        elseif xc == nod(in)
            phi(k,in) = 1;
        end
        k = k + 1;
    end
end

% Half triangles:
k = 1;
for xc = x
    % Right triangle
    if nod(1) < xc && xc < nod(2)
        h2 = nod(2)-nod(1);
        phi(k,1) = (nod(2)-xc)/h2;
    elseif xc == nod(1)
        phi(k,1) = 1;
    end
    % Left triangle
    if nod(nv-1) < xc && xc < nod(nv)
        h1 = nod(nv)-nod(nv-1);
        phi(k,nv) = (xc-nod(nv-1))/h1;
    elseif xc == nod(nv)
        phi(k,nv) = 1;
    end
    k = k + 1;
end

end

