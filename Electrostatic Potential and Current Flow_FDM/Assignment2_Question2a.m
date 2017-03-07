% Assignment 2 - Imran Haider
% Question 2

% Defining Constants
L = 30;
W = 20;
N = 100;              % number of iterations

xvec = linspace(0,L,30);
yvec = linspace(0,W,20);
xvec = xvec';

cMap = zeros(length(xvec),length(yvec));

% conductivity grid cMap
for i = 1:length(xvec)
   for j = 1:length(yvec)
       if (i>10) && (i<20) && (j>0) && (j<8)
           cMap(i,j) = 10^-2;
       elseif (i>10) && (i<20) && (j>12) && (j<21)
           cMap(i,j) = 10^-2;    
       else
           cMap(i,j) = 1;
       end 
   end
end

% Plotting the solution
 [X,Y] = meshgrid(yvec,xvec);
 %surf(X,Y,cMap,'EdgeColor','c');  

nx = length(xvec);
ny = length(yvec);

G = sparse(nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

V = G\B';

vMap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        vMap(i, j) = V(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (vMap(i + 1, j) - vMap(i, j));
        elseif i == nx
            Ex(i, j) = (vMap(i, j) - vMap(i - 1, j));
        else
            Ex(i, j) = (vMap(i + 1, j) - vMap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (vMap(i, j + 1) - vMap(i, j));
        elseif j == ny
            Ey(i, j) = (vMap(i, j) - vMap(i, j - 1));
        else
            Ey(i, j) = (vMap(i, j + 1) - vMap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap .* Ex;
eFlowy = cMap .* Ey;

figure (1)
surf(cMap,'EdgeColor','c');
title('Conductivity Map vs position')
xlabel('Length of region')
ylabel('Width of region')

figure (2)
surf(vMap,'EdgeColor','c');
title('Voltage vs position')
xlabel('Length of region')
ylabel('Width of region')

figure (3)
quiver(Ex', Ey');
title('Electric Field potential vs position')
axis([0 L 0 W]);
xlabel('Length of region')
ylabel('Width of region')

figure (4)
quiver(eFlowx', eFlowy');
title ('Current Density vs position')
axis([0 L 0 W]);
xlabel('Length of region')
ylabel('Width of region')
