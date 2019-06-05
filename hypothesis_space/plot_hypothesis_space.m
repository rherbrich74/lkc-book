% PLOT_HYPOTHESIS_SPACE   plots the hypothesis space with 3 training
%                         examples 
%
% 2001 written by Ralf Herbrich
% Microsoft Research Cambridge
%
% (c) 2001 Microsoft Corporation.  Reproduced with permission.  All rights reserved.

figure;

%% set of training point of unit length 
%% the data points are given by their latitude and azimuth
data = [[0,pi/4]; [0,-pi/4]; [pi/4,pi/4]];

%% plot sphere
[x,y,z] = sphere (50);
s=surfl (x,y,z);
colormap gray;
shading interp;
alpha (0.6);
hold on;

%% plot grand circles
grid = -1.1:0.05:1.1;
for i=1:length (data)
    [gx,gy,gz] = plane (grid, data (i,1), data (i,2));
    h (i) = surfl (gx,gy,gz);
    shading interp;
    [gx,gy,gz] = grandcircle (100, data (i,1), data (i,2), 1.015);
    plot3(gx,gy,gz,'k-','LineWidth', 3);
end

%% one particular weight vector (polar coordinates)
phi = pi*1.35;
theta = 0.1;

%% transform into cartesian coordinates
xs = (cos (phi) * cos (theta)) * 1;
ys = (sin (phi) * cos (theta)) * 1;
zs = (sin (theta)) * 1;

%% plot the additonal weight vector
[xp,yp,zp] = ellipsoid (xs, ys, zs, 0.03, 0.03, 0.03); 
hold on;
p = surfl (xp, yp, zp);
shading interp;

p = plot3 ([0;xs],[0;ys],[0;zs],'b-','LineWidth', 2);
set (p, 'Color', [1, 1, 1]);


% set (p, 'EdgeColor', [0.1 0.1 0.1]);
set (s, 'EdgeColor', [0.3, 0.3, 0.3]);
set (s, 'EdgeAlpha', 0.4);
for i=1:length (data)
%    set (h (i), 'EdgeColor', [0.5, 0.5, 0.5]);
end 

%% finally, tun off the axises and define a nice view
axis off;
axis equal;
vw =[[0.7133 -0.7009  0.0000  -0.0062];
     [0.2049 0.2085   1.3524  -0.8829];
     [0.6703 0.6821   -0.4135  9.5306];
     [0      0        0        1.0000]];
view(vw)

print -deps ../../ps/hypothesis_space.ps
