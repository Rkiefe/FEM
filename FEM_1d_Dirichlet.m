% Author: Rodrigo Kiefe
% Date: 15/02/2023

% Purpose:
% This code simulates the displacement of an elastic bar with Dirichlet
% boundary conditions: Meaning fixed ends

% Recommended Reading:
% Mats G.Larson - Theory Implementation and Applications
% Z. Chen - Finite Element Methods and Their applications

clear
close all
clc

%% Setup
% Bar length
L = 1;          % Length of bar

% Mesh
nv = 50; % number of nodes / vertices

nod = linspace(0,L,nv);

f = 1./(1+20.*nod);

res = 1e-2;
x = 0:res:L;

% Boundary Conditions: u(0) = u(L) = 0

%% Trial function phi:

phi = trial_fun(x,nod);

dphi = zeros(numel(x),nv);
for in = 1:nv
    dphi(:,in) = gradient(phi(:,in))/res;
end

%% Matrix A
A = zeros(nv-2);
k1 = 0; k2 = 0;
for i = 2:nv-1
    k1 = k1+1;
    for j = 2:nv-1
        k2 = k2+1;
        A(k2,k1) = trapz(x,dphi(:,i).*dphi(:,j));
    end
    k2 = 0;
end


%% Right side of equation

F = zeros(nv-2,1);
k = 0;
for i = 2:nv-1
    k = k +1;
    F(k) = trapz(x,f(i)*phi(:,i));
end

%% Solve system of equations:

coef = linsolve(A,F);

% solution:
u = zeros(nv,1);
k = 0;
for in = 2:nv-1
    k = k + 1;
    u(in) = coef(k);
end

%% Plot
subplot(2,1,1)
plot(nod,u,'b.')
xlabel('x'); ylabel('u(x)')
title('1D Elastic Rod')

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

