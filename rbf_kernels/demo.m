%
% DEMO       creates 3D plots of kernel density estimates
%
% 2001 written by Ralf Herbrich
% Microsoft Research Cambridge
%
% (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

figure
density (10, 2, 0.5, 1);
print -deps ps/rbf_density_sigma=0.5.ps
disp ('Press any key for next plot'); pause
density (10, 2, 0.7, 1);
print -deps ps/rbf_density_sigma=0.7.ps
disp ('Press any key for next plot'); pause
density (10, 2, 1.0, 1);
print -deps ps/rbf_density_sigma=1.0.ps
disp ('Press any key for next plot'); pause
density (10, 2, 2.0, 1);
print -deps ps/rbf_density_sigma=2.0.ps
disp ('Press any key for next plot'); pause
