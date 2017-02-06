% Assignment 1 - Question 1
% Imran Haider - 100955365


% Constants
k = 1.380e-23;
W = .10e-6;
L = .2e-6;
mo = 0.26;
T = 300;
m = 9.109e-31*mo;

% Thermal velocity v_th
v_th = sqrt((2*k*T)/m);

% Mean free path
T_mn = 0.2e-12;
d_m = v_th*T_mn;

% Randomized placement of particles
np = 15;
loc = rand(np,2);
xp_loc = loc(:,1)*L;
yp_loc = loc(:,2)*W;

% Plot randomized initial postions
figure (1)
plot(xp_loc,yp_loc,'X')
axis([0 L 0 W])

hold on

% Randomize components of Thermal velocity
vx = 0.5*randn(np,1).*v_th;
vy = 0.5*randn(np,1).*v_th;

mean_vel  = mean(vx.^2 + vy.^2);

dt = 0.05e-13;
nSteps = 500;
for i = 1:nSteps
  
  %t = (m*mean_vel)/(2*k)
    
  % increment x & y positions  
  x_loc = xp_loc  + dt.*vx;  
  y_loc = yp_loc  + dt.*vy;
         
  
  % use logical indexing to find particles approaching the Side-bounds
  ixh = x_loc > L;
  x_loc (ixh) = 0; 
  
  ixl = x_loc < 0;
  x_loc (ixl) = L;
  
  % use logical indexing to find particles approaching Max Y
  iyh = y_loc > W;
  vy(iyh) = -vy(iyh);
  
  % use logical indexing to find particles approaching 0
  iyl = y_loc < 0;
  vy(iyl) = -vy(iyl);
  
  % plot the trajectories
  plot(x_loc,y_loc,'O')
  pause(0.1)
  
  % Update the previous location
  xp_loc = x_loc;
  yp_loc = y_loc;
    
end
%Pixels for Temperature and Density plots
for i = 1:20
    pixelx = W/20; %Dividing the plot into pixels
    
    x1 = (-100e-9)+(i-1)*pixelx;
    x2 = x1 + pixelx;
    
    Ix = x_loc > x1 & x_loc < x2;
    
    for j = 1:20
        pixely = L/20;%Dividing the plot into pixels
        
        y1 = (-50*10^-9)+(j-1)*pixely;
        y2 = y1 + pixely;
        
        Iy = y_loc > y1 & y_loc < y2;
        
        Ixy = Ix & Iy;
        
        xpositionnew = x_loc(Ixy);
        ypa = y_loc(Ixy);
        vxpa = vx(Ixy);
        vypa = vy(Ixy);
        npa = length(xpositionnew);
        
        if npa == 0
            Temperature(i,j) = 300;
        else
            Temperature(i,j) = m*sum(vxpa.^2 + vypa.^2)/(2*k*npa);%rearranged
        end
        
    end
end
    
figure (2)
surf(Temperature)
title('Temperature Plot')
xlabel('Width (m)');
ylabel('Length (m)');
zlabel('Average Temperature (K)');