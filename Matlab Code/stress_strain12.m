function [stress_xy,stress_12,strain_12] = stress_strain12(Qbar_layer,strain_xy,theta)
% Input: Global strain (xy), Qbar, theta
% Output: Global stress, Principal stress, Principal strain
Qb11 = Qbar_layer(1);Qb12 = Qbar_layer(2);Qb16 = Qbar_layer(3);
Qb22 = Qbar_layer(4);Qb26 = Qbar_layer(5);Qb66 = Qbar_layer(6);
Qbar = [Qb11 Qb12 Qb16;...
        Qb12 Qb22 Qb26;...
        Qb16 Qb26 Qb66];
stress_xy = Qbar*strain_xy;
m = cosd(theta);n = sind(theta);
T = [m^2 n^2 2*m*n;...
     n^2 m^2 -2*m*n;...
     -m*n m*n m^2-n^2];
stress_12 = T*stress_xy;
strain_12 = T*[strain_xy(1);strain_xy(2);strain_xy(3)/2];
strain_12(3) = 2*strain_12(3);
end

