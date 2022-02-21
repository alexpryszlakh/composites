% Lamina Stiffness Analyzer
% Author: Daniel Donate-Perez
% Date:   3/8/2021
clear
close all
clc
format loose
format shortG
%% Effective Elastic Constants
% Inputs
% E1 = input('Elastic Modulus in Fiber Direction, E1 [Pa]: ');
% E2 = input('Elastic Modulus in Transverse Direction, E2 [Pa]: ');
% G12 = input('Shear Modulus in 1-2 Plane, G12 [Pa]: ');
% v12 = input('Poisson Ratio in 1-2 Plane, v12: ');
% fprintf('\n')
% fprintf('Enter Lamina Orientation, in Global Coordinate System x-y-z, CCW+ \n \r')
% theta = input('Lamina Orientation, theta [degrees]: ');%orientation [degrees]
% Calculation
% s = sind(theta);c = cosd(theta); %
% Ex = 1/( c^4/E1 + (1/G12 - 2*v12/E1)*s^2*c^2 + s^4/E2 ); %
% Ey = 1/( s^4/E1 + (1/G12 - 2*v12/E1)*s^2*c^2 + c^4/E2 ); %
% Gxy = 1/( (s^4 + c^4)/G12 + 4*(1/E1 + 1/E2 + 2*v12/E1 - 1/(2*G12))*s^2*c^2 ); %
% vxy = Ex*( v12*(s^4 + c^4)/E1 - (1/E1 + 1/E2 - 1/G12)*s^2*c^2 ); %
%% Loop
GF_prop = [6.53 2.6668 0.945520205839657 0.254];
CF_prop = [22.7 1.1328 0.5061 0.2675]; 
theta_vect = 0:1:360;
Stiff_Matrix_CF = zeros(length(theta_vect),5);
Stiff_Matrix_GF = zeros(length(theta_vect),5);
for i = 1:length(theta_vect)
    CF_prop(5) = theta_vect(i);
    GF_prop(5) = theta_vect(i);
    Stiff_Matrix_CF(i,1:5) = LamEffElastic(CF_prop);
    Stiff_Matrix_GF(i,1:5) = LamEffElastic(GF_prop);
end
subplot(2,2,1)
plot(theta_vect,Stiff_Matrix_CF(:,1))
hold on
plot(theta_vect,Stiff_Matrix_GF(:,1))
title('Ex')
subplot(2,2,2)
plot(theta_vect,Stiff_Matrix_CF(:,2))
hold on
plot(theta_vect,Stiff_Matrix_GF(:,2))
title('Ey')
subplot(2,2,3)
plot(theta_vect,Stiff_Matrix_CF(:,3))
hold on
plot(theta_vect,Stiff_Matrix_GF(:,3))
title('Gxy')
figure(2)
subplot(1,2,1)
polarplot(theta_vect*((2*pi)/360),Stiff_Matrix_GF(:,1),...
          theta_vect*((2*pi)/360),Stiff_Matrix_CF(:,1))
title('E_1 [msi]')
legend('Glass Fiber Reinforced','Carbon Fiber Reinforced','Location','southoutside')
pax = gca;
thetaticks(0:45:90);

subplot(1,2,2)
polarplot(theta_vect*((2*pi)/360),Stiff_Matrix_GF(:,3),...
          theta_vect*((2*pi)/360),Stiff_Matrix_CF(:,3))
          
thetaticks(0:45:90);
title('G_1_2 [msi]')
legend('Glass Fiber Reinforced','Carbon Fiber Reinforced','Location','southoutside')
% subplot(2,2,4)
% plot(theta_vect,Stiff_Matrix_CF(:,4))
% hold on
% plot(theta_vect,Stiff_Matrix_GF(:,4))
% title('vxy')    
