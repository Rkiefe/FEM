function itg = bordintg(cb,nd1,nd2,p,t)

[x1,x2,y1,y2,m1,m2,phi1,phi2,cbnd] = egde(nd1,cb,p,t);
[x3,x4,y3,y4,m3,m4,phi3,phi4,cbnd2] = egde(nd2,cb,p,t);

f = [phi1,phi2]; f = f(not(isnan(f)));
g = [phi3,phi4]; g = g(not(isnan(g)));
xc1 = [x1,x2]; xc1 = xc1(not(isnan(xc1)));
xc2 = [x3,x4]; xc2 = xc2(not(isnan(xc2)));
yc1 = [y1,y2]; yc1 = yc1(not(isnan(yc1)));
yc2 = [y3,y4]; yc2 = yc2(not(isnan(yc2)));

% \\ ---------- Intersection of edges -------- //

[xt,i1,i2] = intersect(xc1,xc2);
[yt,j1,j2] = intersect(yc1,yc2);


for ind = 1:2
    cx = eval(sprintf('x%1d',ind));
    cy = eval(sprintf('y%1d',ind));
    if isequal(sort(xt),sort(cx))
        m = eval(sprintf('m%1d',ind));
        i = ind;
    end
    if isequal(sort(yt),sort(cy))
        m = eval(sprintf('m%1d',ind));
        i = ind;
    end
end
for ind = 3:4
    cx = eval(sprintf('x%1d',ind));
    cy = eval(sprintf('y%1d',ind));
    if isequal(sort(xt),sort(cx))
        m = eval(sprintf('m%1d',ind));
        j = ind;
    end
    if isequal(sort(yt),sort(cy))
        m = eval(sprintf('m%1d',ind));
        j = ind;
    end
end

% \\ ---------- Integral -------- //

% if nd1 and nd2 are different nodes
if nd1~=nd2
    % if they are well behaved nodes
    if numel(xt) >= 1 && numel(yt) >= 1
        % But such that they connect only on one point
        if numel(xt) == 1 && numel(yt) == 1
            itg = 0;
            return
        end

        % else, trully well behaved
        f = eval(sprintf('phi%1d',i));
        xf = eval(sprintf('x%1d',i));
        yf = eval(sprintf('y%1d',i));

        g = eval(sprintf('phi%1d',j));
        xg = eval(sprintf('x%1d',j));
        yg = eval(sprintf('y%1d',j));

        if isequal(f,g)
            g = flip(g);
        end

        if isnan(m)
            itg = trapz(yt,f.*g);
        else
            itg = trapz(xt,f.*g)*(1+m^2);
        end

        % if they are space too far appart ...
    else
        % ... such that no edges connect
        itg = 0;
    end

    % If both nodes are the same node
else
    if isnan(phi2)
        if isnan(m1)
            itg = trapz(y1,phi1.*phi1);
        else
            itg = trapz(x1,phi1.*phi1).*sqrt(1+m1^2);
        end
    elseif isnan(m1) && isnan(m2)
        itg = trapz(y1,phi1.*phi1) + trapz(y2,phi2.*phi2);
    elseif isnan(m1)
        itg = trapz(y1,phi1.*phi1) + trapz(x2,phi2.*phi2).*sqrt(1+m2^2);
    elseif isnan(m2)
        itg = trapz(x1,phi1.*phi1).*sqrt(1+m1^2) + trapz(y2,phi2.*phi2);
    else
        itg = trapz(x1,phi1.*phi1).*sqrt(1+m1^2) + trapz(x2,phi2.*phi2).*sqrt(1+m2^2);
    end
end

end
