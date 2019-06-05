% PLOT_DATA_SPACE   plots the data space with 3 training examples 
%
% 2001 written by Ralf Herbrich
% Microsoft Research Cambridge
%
% (c) 2001 Microsoft Corporation.  Reproduced with permission.  All rights reserved.

figure;

%% set of training point of unit length 
%% the data points are given by their latitude and azimuth
data = [[0,pi/4]; [0,-pi/4]; [pi/4,pi/4]];

%% one particular weight vector (polar coordinates)
phi = pi*1.35;
theta = 0.1;

% plot the data points by little blobs;
for i=1:length (data)
  xp = 1.4 * (cos (data (i,1)) * cos (data(i,2)));
  yp = 1.4 * (sin (data (i,1)) * cos (data(i,2)));
  zp = 1.4 * (sin (data(i,2)));

  [x,y,z] = ellipsoid (xp, yp, zp, 0.05, 0.05, 0.05);
  surfl (x,y,z);
  shading interp;
  hold on;
  
  set (plot3 ([0; xp], [0; yp], [0;zp], 'w-', 'LineWidth', 2), 'Color', [0.9 0.9 0.9]);
end

%% plot the hyperplane 
hold on;
[gx,gy,gz] = plane (-2:0.2:2, phi, theta);
s = surfl (gx,gy,gz);
shading interp;
alpha (0.4);
colormap gray;

%% create meshgrids
set (s, 'EdgeColor', [0.5, 0.5, 0.5]);
set (s, 'EdgeAlpha', 0.6);

%% draw the coordinate axes 
set (plot3 ([-2.8;+2.8],[0;0],[0;0], 'k-', 'LineWidth', 1), 'Color', [0.2 0.2 0.2]);
set (plot3 ([0;0],[-2.3;+2.3],[0;0], 'k-', 'LineWidth', 1), 'Color', [0.2 0.2 0.2]);
set (plot3 ([0;0],[0;0],[-2.3;+2.3], 'k-', 'LineWidth', 1), 'Color', [0.2 0.2 0.2]);

%% finally, tun off the axises and define a nice view
axis off;
axis equal;

VW = [[0.5015   -0.8652    0.0000    0.1819];
      [0.2443    0.1416    0.9593   -0.6726];
      [0.8300    0.4811   -0.2824    8.1459];
      [     0         0         0    1.0000]];
view (VW);
zoom (1.5);

print -deps ../../ps/data_space.ps
