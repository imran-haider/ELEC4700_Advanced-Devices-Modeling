% Assignment 2 - Imran Haider
% Question 1

% Define constants
L = 30;
W = 20;
N = 20;              % number of iterations

% Defining x and y vectors
numx = 60;
% dx = 2/(numx-1);
x = linspace(0,L,N);

numy = 60;
dy = 2/(numy-1);
y = 0:1:W;

% Defining placeholder Zero matrices
v = zeros(length(x),length(y));
% vn = zeros(numx,numy);

%Boundary conditions
%v(:,1)=0;
%v(:,numx)=y;
v(1,:)= 5;
v(end,:)=0;

% j=2:numx-1;
% i=2:numy-1;

for in = 1 : N
    for j=1:length(y)
        for i=2:length(x)-1
            if (j==21)
                v(i,j)=(v(i,j-1));%No V_mat(idx,ydx+1) since there is no dot above it
            elseif(j==1)
                v(i,j)=(v(i,j+1));%No +V_mat(idx,ydx-1) since there is no dot bellow it
            else
                v(i,j)=(v(i+1,j)+v(i-1,j)+v(i,j-1)+v(i,j+1))/4;
            end
        end
    end
end

% Plotting the solution
[Y,X] = meshgrid(y,x);
surf(X,Y,v,'EdgeColor','c');       
colorbar
title('Finite Difference Method with fixed Boundary Conditions')
xlabel('Length of the region')
ylabel('Width of the region')
zlabel('Voltage (V)')

%     vn=v;
%     v(i,j)=((dy^2*(vn(i+1,j)+vn(i-1,j)))+(dx^2*(vn(i,j+1)+vn(i,j-1))))/(2*(dx^2+dy^2));
%     v(:,1)=1;
%     v(:,numx)=0;
% %     %v(1,:)=0;
% %     %v(numy,:)=0;   


