clear
close all
clc
load('grid.mat')
load('trial.mat')

xs = linspace(0,100,200);
ys = linspace(0,50,200);
vert = [45,22.5;55,22.5;55,27.5;45,27.5];

nb = labels(end);

%% Borders and conditions

% Nod indices
nd_ind = 1:nv;

% Conditions
jcurr = zeros(numel(xs),numel(ys));

mu0 = pi*4e-7;

mu  = mu0.*ones(numel(xs),numel(ys));
[in,on] = inpolygon(xs,ys,vert(:,1),vert(:,2));  % metal
mu(in,in) = 5000.*(pi*4e-7);

% border indices at which g and gam are ~= 0
gb = [1,3,5,6,7,8];
gamb = gb;
g = [-10,10,0,0,0,0];   % value of g on those borders
gam = [1, 1,0,0,0,0];    % value of gam on those borders

%% Gradient of phi
dpsi_dx = cell(nv,1);
dpsi_dy = cell(nv,1);

for i = 1:nv
    [dpsi_dy{i},dpsi_dx{i}] = gradient(psi{i});
end

%% (M+R)c = b+s

M = zeros(nv);
R = M;
b = zeros(nv,1);
s = b;

% M
for in1 = 1:nv
    nd1 = nd_ind(in1);
    for in2= 1:nv
        nd2 = nd_ind(in2);
        f = (dpsi_dx{nd1}.*dpsi_dx{nd2} + dpsi_dy{nd1}.*dpsi_dy{nd2})./mu;
        M(in1,in2) = trapz(ys,trapz(xs,f));
    end
end

% b
for in1 = 1:nv
    nd = nd_ind(in1);
    f = psi{nd}.*jcurr;
    b(in1) = trapz(ys,trapz(xs,f));
end

% R
%      f = gam.*psi{nd1}.*psi{nd2};
for i1 = 1:numel(bnod)
    nd1 = bnod(i1);
    cb1 = p(3,nd1);
    i = find(nd_ind==nd1);
    for i2 = 1:numel(bnod)
        nd2 = bnod(i1);
        cb2 = p(3,nd1);
        j = find(nd_ind==nd2);

        [~,chk] = intersect(gamb,cb1); % if cb1 is in gamb, then the itg ~= 0

        if cb1 == cb2 && not(isempty(chk))
            itg = bordintg(cb1,nd1,nd2,p,t);
            R(i,j) = abs(itg)*gam(chk)/mu0;
        end
    end
end


% s
%     f = psi{nd}.*g;
for ind = 1:numel(bnod)
    nd = bnod(ind);
    cb = p(3,nd);
    j = find(nd_ind==nd);

    [~,chk] = intersect(gb,cb); % if cb1 is in gamb, then the itg ~= 0

    if not(isempty(chk))
        [x1,x2,y1,y2,m1,m2,phi1,phi2,cbnd] = egde(nd,cb,p,t);

        if not(isnan(m1)) && not(isnan(m2)) % regular well behaved case
            itg = trapz(x1,phi1).*sqrt(1+m1^2) + trapz(x2,phi2).*sqrt(1+m2^2);

        else
            if isnan(phi2)  % single edge case
                if isnan(m1)
                    itg = trapz(y1,phi1);
                else
                    itg = trapz(x1,phi1).*sqrt(1+m1^2);
                end

                % if not a single edge ...
            elseif isnan(m1) && not(isnan(m2))
                itg = trapz(y1,phi1) + trapz(x2,phi2).*sqrt(1+m2^2);
            elseif isnan(m2) && not(isnan(m1))
                itg = trapz(y2,phi2) + trapz(x1,phi1).*sqrt(1+m1^2);
            else
                itg = trapz(y2,phi2) + trapz(y1,phi1);
            end
        end
        s(j) = abs(itg)*g(chk)/mu0;
    end
end


%% Linsolve

coef = linsolve(M+R,s+b);

Ap = zeros(numel(xs),numel(ys));
for ind = 1:nv
    Ap = Ap + coef(ind)*psi{nd_ind(ind)};
end

%% Plot vector Potential A
figure
surf(xs,ys,Ap',"EdgeColor","none")
title('Vector Potential')

%% Magnetic Field Lines
% B = rot A -> dA/dy i - dA/dx j

[dAy,dAx] = gradient(Ap);

Bx = dAy;
By = -1.*dAx;

%% Plot Magnetic Field
figure
quiver(xs,ys,Bx',By'); hold on

% Add borders
for ind = 1:nb
    cb = eval(sprintf('b%1d',ind)); % current border
    plot(cb(:,1),cb(:,2),'r-')
end
xlabel('x'); ylabel('y')
title('Magnetic Field')

%% Plot Magnetic field intensity
figure
Ix = Bx'; Iy = By';
surf(xs,ys,sqrt(Ix.^2+Iy.^2),'edgecolor','none')
title('Magnetic field Strength (T)')
view(2)
colorbar
%% Save
% save('Result.mat')

%% Functions
function plotri(p,t,i,style)
if nargin < 4
    style = 'b.-';
end
plot([p(1,t(1,i)),p(1,t(2,i)),p(1,t(3,i)),p(1,t(1,i))],[p(2,t(1,i)),p(2,t(2,i)),p(2,t(3,i)),p(2,t(1,i))],style)
end

function pf = findnode(p,px,py)

d = 2000; % arbitrarly big
for j = 1:length(p)

    d_n = sqrt((p(1,j)-px)^2 + (p(2,j)-py)^2 );

    if d_n < d
        d = d_n;
        pf = j;
    end
end

end

function t_ind = find_triang(pf,t)
t_ind = [];

for j = 1:length(t)
    for k = 1:3
        if pf == t(k,j)
            t_ind(end+1) = j;
        end
    end
end

end

function logic = isequiv(a,b)
logic = true;
if sum(a) ~= sum(b)
    logic = false;
    return
else
    for i = 1:length(a)
        if a(i)~= b(i)
            logic = false;
            return
        end
    end
end
end

function area = t_area(t,p)

for i = 1:length(t)
    pt(:,i) = p(:,t(i));
end

x = pt(1,:); y = pt(2,:);
area = polyarea(x,y);
end

function ib = yind(ss,b)
ib = [];
for i1 = 1:numel(b)
    d = 999;
    for i2 = 1:numel(ss)
        if abs(ss(i2)-b(i1)) < d
            aux = i2;
            d = abs(ss(i2)-b(i1));
        end
    end
    ib(numel(ib)+1) = aux;
end
end




