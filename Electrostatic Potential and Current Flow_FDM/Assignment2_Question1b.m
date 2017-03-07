% Assignment 2 - Imran Haider
% Question 1

% Define constants
close all; clc; clear;
a = 2;
b = 1;
v0 = 1;
N = 100;              % number of iterations

numx = 30;
dx = 2/(numx-1);
x = -b:dx:b;

numy = 60;
dy = 2/(numy-1);
y1 = 0:dy:a;
y=y1';

v = zeros(length(y),length(x));
%vnew = zeros(length(x),length(x));

c = (4*v0)/pi;

for n = 1:2:N
    vnew = c*(1/n).*((cosh(n*pi.*x./a)./cosh(n*pi*b/a)).*(sin(n*pi.*y./a)));
     v = v + vnew;
     v(:,1)=v0;
     v(:,end)=v0;
     v(1,:)=0;
     v(end,:)=0;
end

% Plotting the solution
surf(x,y,v,'EdgeColor','c');       
colorbar
title({'Analytical Solution';['{\itNumber of iterations} = ',num2str(N)]})
xlabel('x-axis')
ylabel('y-axis')
zlabel('Solution (v)') 


