% Assignment 1 - Question 3
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
np = 10;
valid = false;
while (valid ~= true)
    loc = rand(np,2);
    xp_loc = loc(:,1)*L;
    yp_loc = loc(:,2)*W;
    x_bounds = xp_loc >= 0.08e-6;
    if (xp_loc(x_bounds) <= (0.04e-6 + 0.08e-6))
        if (yp_loc(x_bounds)> 0.04e-6) 
            if (yp_loc(x_bounds) < 0.06e-6)
            valid = true;
            end
        end
    else
        valid = false;
    end
end

% Plot randomized initial postions
figure (1)
plot(xp_loc,yp_loc,'X')
axis([0 L 0 W])

hold on

% Randomize components of Thermal velocity
vx = 0.5*randn(np,1).*v_th;
vy = 0.5*randn(np,1).*v_th;

dt = 0.5e-14;
nSteps = 500;
mVel = zeros(nSteps,1);
pscat = 1 - exp(-dt/T_mn);

% Making the boxes
pos_bottom = [0.08e-6 0  0.04e-6 0.04e-6];
rec = rectangle('Position',pos_bottom)
rec.LineWidth = 0.1;

pos_top = [0.08e-6 0.06e-6  0.04e-6 0.04e-6];
rec = rectangle('Position',pos_top)
rec.LineWidth = 0.1;

% Upper Box Bounds
topBoxLeftSide = 0.08e-6;
topBoxRightSide = (0.08e-6 + 0.04e-6);
topBoxUp = W;
topBoxDown = 0.06e-6;

% Bottom Box Bounds
lowBoxLeftSide = 0.08e-6;
lowBoxRightSide = (0.08e-6 + 0.04e-6);
lowBoxUp = 0.04e-6;
lowBoxDown = 0;

for i = 1:nSteps
    
  %Scattering
  r = rand([np,1]);   
  p = pscat > r;
    
  vx(p) = randn().*v_th;
  vy(p) = randn().*v_th;
    
  % increment x & y positions  
  x_loc = xp_loc  + dt.*vx;  
  y_loc = yp_loc  + dt.*vy;
         
  
  % use logical indexing to find particles approaching the Side-bounds
  ixh = x_loc > L;
  x_loc (ixh) = x_loc(ixh) - L;
  
  ixl = x_loc < 0;
  x_loc (ixl) = x_loc(ixl) + L;
  
  % use logical indexing to find particles approaching the sides of Boxes 
  ixLeftSide = x_loc > lowBoxLeftSide & xp_loc < lowBoxLeftSide & y_loc < lowBoxUp;
  x_loc(ixLeftSide) = x_loc(ixLeftSide) - ((x_loc(ixLeftSide) - lowBoxLeftSide)*2);
  vx(ixLeftSide) = -vx(ixLeftSide);
  
  ixRightSide = x_loc < lowBoxRightSide & xp_loc > lowBoxRightSide & y_loc < lowBoxUp;
  x_loc(ixRightSide) = ((lowBoxRightSide - x_loc(ixRightSide))*2) + x_loc(ixRightSide);
  vx(ixRightSide) = -vx(ixRightSide);
      
  % use logical indexing to find particles approaching the sides of Boxes 
  ixTopLeftSide = x_loc > topBoxLeftSide & xp_loc < topBoxLeftSide & y_loc > topBoxDown;
  x_loc(ixTopLeftSide) = x_loc(ixTopLeftSide) - ((x_loc(ixTopLeftSide) - topBoxLeftSide)*2);
  vx(ixTopLeftSide) = -vx(ixTopLeftSide);

  
  ixTopRightSide = x_loc < topBoxRightSide & xp_loc > topBoxRightSide & y_loc > topBoxDown;
  x_loc(ixTopRightSide) = ((topBoxRightSide - x_loc(ixTopRightSide))*2) + x_loc(ixTopRightSide);
  vx(ixTopRightSide) = -vx(ixTopRightSide);

  
  % use logical indexing to find particles approaching Max Y
  iyh = y_loc > W;
  y_loc(iyh) = y_loc(iyh) - ((y_loc(iyh) - W)*2);
  vy(iyh) = -vy(iyh);
  
  % use logical indexing to find particles approaching 0
  iyl = y_loc < 0;
  y_loc(iyl) = ((0 - y_loc(iyl))*2) + y_loc(iyl);
  vy(iyl) = -vy(iyl);
  
  % checking the center area for boundary conditions
  lbt = x_loc < lowBoxRightSide & x_loc > lowBoxLeftSide & y_loc <= lowBoxUp;
  y_loc(lbt) = ((lowBoxUp - y_loc(lbt))*2) + y_loc(lbt);
  vy(lbt) = -vy(lbt);
  
  tbd = x_loc < lowBoxRightSide & x_loc > lowBoxLeftSide & y_loc >= topBoxDown;
  y_loc(tbd) = y_loc(tbd)-(2*(y_loc(tbd)-topBoxDown));
  vy(tbd) = -vy(tbd);
  
  % plot the trajectories
    plot(x_loc,y_loc,'O','MarkerSize',2)
    title('Random Motion with Scattering');
    xlabel('Length(m)');
    ylabel('Width(m)');
    pause(0.1)
    grid on
    
  % Update the previous location
  xp_loc = x_loc;
  yp_loc = y_loc;
    
end
%Pixels for Temperature and Density plots
for i = 1:20
    pixelx = W/20; %Dividing the plot into pixels
    
    x1 = (-100*10^-9)+(i-1)*pixelx;
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
            elecdensity(i,j) = npa./(W*L);
        else
            Temperature(i,j) = m*sum(vxpa.^2 + vypa.^2)/(2*k*npa);%rearranged
            %equation for KT=(mv^2)/2
            elecdensity(i,j) = npa./(W*L);
        end
        
    end
end
    
figure(2)
surf(elecdensity)
title('Electricon Density map')
xlabel('Width (m)');
ylabel('Length (m)');
zlabel('Electron Density per Pixal');

figure (3)
surf(Temperature)
title('Temperature Plot')
xlabel('Width (m)');
ylabel('Length (m)');
zlabel('Average Temperature (K)');