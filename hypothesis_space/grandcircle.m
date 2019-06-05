% GRANDCIRCLE    computes the coordinates of the grandcircle given
%                normal vector (phi,theta)
%
% USAGE: [x,y,z] = grandcircle (N,phi,theta)
% 
% parameters: N         number of points of the grand circle
%             phi       rotation of plane around z axis
%             theta     rotation of plane around x axis
%             x,y,z     coordinates of the grandcircle
%
% (c) 2001 Microsoft Corporation.  Reproduced with permission.  All rights reserved.

function [x,y,z] = grandcircle (N,phi,theta,off)

  x = zeros (N);
  y = zeros (N);
  z = zeros (N);
  
  for i=1:N
    alpha = 2*pi/(N-1)*(i-1);
    beta = -atan2 (cos (theta), sin (phi) * sin (theta) * sin (alpha) + cos ...
		   (phi) * sin (theta) * cos (alpha));
    
    x (i) = sin (alpha) * sin (beta) * off;
    y (i) = cos (alpha) * sin (beta) * off;
    z (i) = cos (beta) * off;
  end

