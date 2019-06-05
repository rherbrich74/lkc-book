% PLANE    generates a hyperplane for a given meshgrid 
%
% USAGE: [x, y, z] = plane (grid, phi, theta)
% 
% parameters: grid      vector of grid points
%             phi       rotation of plane around z axis
%             theta     rotation of plane around x axis
%             x         x-coordinates of the mesh grid
%             y         y-coordinates of the mesh grid
%             z         z-coordinates of the mesh grid
%
% (c) 2001 Microsoft Corporation.  Reproduced with permission.  All rights reserved.

function [x, y, z] = plane (grid, phi, theta)
  
  N = length (grid);
  
  for i=1:N
    for j=1:N
      x (i,j) = grid (i);
      y (i,j) = grid (j);
      z (i,j) = -sin (phi) * tan (theta) * x (i,j) - cos (phi) * tan (theta) ...
		* y(i,j); 
    end
  end
  
