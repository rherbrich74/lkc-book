% DENSITY   plots a desnity estimate of the input distribution
%           using RBF kernels
%
%      DENSITY(N,M) makes a denisty plot of N data-points with a 
%         distance of M and an RBF kernel with SIGMA=1
%
%      DENSITY(N,M,SIGMA) makes a denisty plot of N data-points with a 
%         distance of M and an RBF kernel with SIGMA
%
%      DENSITY(N,M,SIGMA,alphaScaling) makes a denisty plot of N data-points with a 
%         distance of M and an RBF kernel with SIGMA=1; the transparency
%         is increased for alphaScaling > 1 and decreased for alphaScaling < 1
%
% 2001 written by Ralf Herbrich
% Microsoft Research Cambridge
%
% (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

function density (N,M,SIGMA,alphaScaling)

    %% parameter check
    if nargin < 3
        SIGMA = 1.0;
    end
    
    if nargin < 4
        alphaScaling = 1.0
    end

    hold off;

    %% generate data (or load, if already generated
    if (exist ('dens.dat') == 2)
        load -ascii 'dens.dat';
        data_x = dens (:,1);
        data_y = dens (:,2);
        clear dens;
    else    
        data_x = [randn(N, 1) - M, randn(N ,1) + M];
        data_y = [randn(N, 1) + M, randn(N ,1) - M];
        save -ascii 'dens.dat' data_x data_y;
    end
    

    rg_x = [min(min(data_x*1.5)) max(max(data_x*1.5))];
    rg_y = [min(min(data_y*1.5)) max(max(data_y*1.5))];
    gr_x = rg_x(1):(rg_x(2)-rg_x(1))/100:rg_x(2);
    gr_y = rg_y(1):(rg_y(2)-rg_y(1))/100:rg_y(2);

    %% generate the surfaces for each single point
    [x,y] = meshgrid(gr_x, gr_y);
    z = zeros (size (x));
    for i=1:2*N
        z = z + exp ((-(x - data_x (i)).^2 - (y - data_y (i)).^2) ./ (2*SIGMA^2));
    end 

    %% plot the surface
    s = surfl (x,y,z);
    shading interp;
    alpha = 1 - exp (-z/alphaScaling);
    colormap gray;
    hold on;
    contour(x,y,z,10);

    %% plot the data points as well
    set (plot (data_x,data_y, 'k.', 'MarkerSize', 12), 'Color', [0.2 0.2 0.2]);
    end
        set (s, 'FaceAlpha', 'interp', 'AlphaDataMapping', 'scale', 'AlphaData', alpha);

    axis normal
    axis tight
