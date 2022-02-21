% Effective Lamina Properties
% Author: Daniel Donate-Perez
% Date:   3/8/2021
clear all, clc
format loose
format shortG
%% Elastic Constants
% Inputs
units = input('PSI (1) or Pa (2) (enter 1 or 2): ');
Vf = input('Fiber Volume Fraction, Vf: ');
fprintf('\n')
fprintf('Enter Fiber Properties \r \n')
E1f = input('Elastic Modulus in Fiber Direction, E1f [Pa]: ');
E2f = input('Elastic Modulus in Transverse Direction, E2f [Pa]: ');
G12f = input('Shear Modulus in 1-2 Plane, G12f [Pa]: ');
v12f = input('Poisson Ratio in 1-2 Plane, v12f: ');
fprintf('\n')
fprintf('Enter Matrix Properties \r \n')
Em = input('Elastic Modulus of Matrix, Em [Pa]: ');
vm = input('Poisson Ratio of Matrix, vm: ');
Gm = Em/(2*(1 + vm));
% % Calculations
E1 = E1f*Vf + Em*(1-Vf); % Rule of Mixtures
v12 = v12f*Vf + vm*(1-Vf); % Rule of Mixtures
E2 = 1/( (1 - sqrt(Vf))/Em + sqrt(Vf)/...
         (E2f*sqrt(Vf) + Em*(1 - sqrt(Vf))) );% Improved Mechanics of Materials model
G12 = Gm*( ((Gm + G12f) - Vf*(Gm - G12f))/...
           ((Gm + G12f) + Vf*(Gm - G12f)) );% CCM
% Outputs
lamina_properties = [E1f,E2f,G12f,v12f;Em,Em,Gm,vm;E1,E2,G12,v12];
if units == 1
    lamina_properties(:,1:3) = lamina_properties(:,1:3)*(1e-6);
elseif units == 2
    lamina_properties(:,1:3) = lamina_properties(:,1:3)*(1e-9);
end
lamina_properties = array2table(lamina_properties,'VariableNames',...
                         {'E1 [msi]','E2 [msi]','G12 [msi]','v12'},...
                          'RowNames',{'Fiber','Matrix','Composite Lamina'})
%% Strength Properties
