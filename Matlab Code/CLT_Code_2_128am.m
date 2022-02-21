% CLT Code 2.0
% Author(s): Daniel Donate-Perez
% Date: 3/17/2021
clear
close all
clc
format shortG
%% Part 0: Material Properties (Done)
% Input:
%   1. Layer material properties
%       a. E1,E2,v12,G12
%       *b. CTE(1,2,3),CMD(1,2,3)*
%   2. Layer geometry
%       a. Number of layers
%       b. Layer thickness
%       c. Layer orientation
% Pre-Output:
%   1. z-vector
%   2. *Record whether symmetric*
%   3. *Record whether balanced*
%   4. *Record whether cross-ply*
% Output:
%   1. Table of layer information 

% Input_1.a: Layer material properties(E1,E2,v12,G12)
CFRP_mech = [155e9 12.1e9 0.248 4.4e9];
% Input_2.a: number of layers
number_layers = input('Number of Layers: ');
% Input_2.b: layer thickness
t = 0.15e-3;
% Input_2.c: layer orientation
count = 1;
layer_orientation = zeros(1,number_layers);
while count <= number_layers
    fprintf('Orientation of layer %.0f: \n',count);
    layer_orientation(count) = input('Layer Orientation [degrees]: ');
    count = count + 1;
end
% Pre-Output.1: z-vector
count = 1;
z_vect = zeros(1,number_layers+1);
z_vect(1) = -t*number_layers/2;
while count <= number_layers
        z_vect(count + 1) = z_vect(1) + t*count;
        count = count + 1;
end
% Pre-Output.2: symmetric laminate
symmetric = input('Is this laminate symmetric? (1=Yes, 0=No): ');
% Pre-Output.3: balanced laminate
balanced = input('Is this laminate balanced? (1=Yes, 0=No): ');
% Pre-Output.4: cross-ply laminate
crossply = input('Is this laminate crossply? (1=Yes, 0=No): ');
% Output: Table of layer information
colNames = {'Layer_Number','E1','E2','v12','G12',...
            'Layer_Orientation','Layer_Thickness'};
laminate(1:number_layers,1) = 1:number_layers;
laminate(:,7) = t;
for i = 1:number_layers
    laminate(i,6) = layer_orientation(i);
    laminate(i,2:5) = CFRP_mech;
end
laminate_table = array2table((laminate),'VariableNames',colNames);    
%% Part 1: Calculate Force-Moment Resultants (Done)
% Input:
%   1. Applied loads
%       a. P(x,y),V(xy),M(x,y,xy)
%   2. Laminate geometry
%       a. Laminate width (along y-axis)
%       b. Laminate length (along x-axis)
% Output:
%   1. Force-Moment Resultants

% Input.1: Applied loads
FM(1) = input('Applied normal force in x-direction: ');
FM(2) = input('Applied normal force in y-direction: ');
FM(3) = input('Applied shear force in xy-plane: ');
FM(4) = input('Applied moment about y-axis: ');
FM(5) = input('Applied moment about x-axis: ');
FM(6) = input('Applied twisting moment about xy-plane: ');
% Input.2: Laminate geometry
y = input('Plate width (along y-axis): ');
x = input('Plate length (along x-axis): ');
% Output: Force-Moment Resultants
NM(1) = FM(1)/y;NM(2) = FM(2)/x;NM(3) = FM(3)/y;
NM(4) = FM(4)/y;NM(5) = -FM(5)/x;NM(6) = -FM(6)/y;
NM = NM';
%% Part 2: Calculate ABD Matrix (Done)
% Input:
%   1. Table of layer information
%   2. z-vector
% Pre-output:
%   1. Table with Qbar
% Output:
%   1. ABD Matrix
%   2. abd Matrix

% Pre-output: Qbar
Qb(1:number_layers,1) = 1:number_layers;
for i = 1:number_layers
    E1 = laminate_table.E1(i);E2 = laminate_table.E2(i);
    v12 = laminate_table.v12(i);G12 = laminate_table.G12(i);
    theta = laminate_table.Layer_Orientation(i);
    Qb(i,2:7) = Qbar(E1,E2,v12,G12,theta);
end

% Output.1: ABD Matrix
ABDmat = ABD(Qb(:,2:7),z_vect);
% If symmetric, Bij = 0
if symmetric == 1
    ABDmat(1:3,4:6) = 0;
    ABDmat(4:6,1:3) = 0;
end
% If balanced, A16 = A26 = 0
if balanced == 1
    ABDmat(1:2,3) = 0;
    ABDmat(3,1:2) = 0;
end
% If cross-ply, A16 = A26 = B16 = B26 = D16 = D26 = 0
if crossply == 1
    ABDmat(1:2,3) = 0; %A16,A26
    ABDmat(3,1:2) = 0; % 
    ABDmat(1:2,6) = 0; %B16,B26
    ABDmat(6,1:2) = 0; %
    ABDmat(4:5,3) = 0; %
    ABDmat(3,4:5) = 0; %
    ABDmat(4:5,6) = 0; %D16,D26
    ABDmat(6,4:5) = 0; %
end

% Output.2: abd Matrix
% abd = inv(ABDmat);
%% Part 3: Calculate Layer Strains
% Input:
%   1. Force-moment resultants
%   2. abd matrix
%   3. Table with layer S,Qbar,T
% Pre-output:
%   Neutral plane strains and Neutral plane curvatures
% Output
%   1. Table with layer strains (x,y,xy)

% Pre-output:
e_k = ABDmat\NM;
strain = zeros(2*number_layers,7);
stress = zeros(2*number_layers,7);
zvect_layer = zeros(2 + 2*(length(z_vect)-2),1);
zvect_layer(1) = z_vect(1);
zvect_layer(end) = z_vect(end);
for i = 3:length(zvect_layer)-1
    index = fix(i/2 + 0.5);
    zvect_layer(i-1:i) = z_vect(index);
end
for i = 1:2*number_layers
    index = fix(i/2 + 0.5);
    strain(i,1) = index;
    ex = e_k(1) + zvect_layer(i)*e_k(4);
    ey = e_k(2) + zvect_layer(i)*e_k(5);
    exy = e_k(3) + zvect_layer(i)*e_k(6);
    strain(i,2:4) = [ex ey exy];
    evect = [ex ey exy]';
    theta = laminate_table.Layer_Orientation(index);
    Qb1 = Qb(:,2:7);
    [stress_xy,stress_12,strain_12] = stress_strain12(Qb1(index,:),evect,theta);
    strain(i,5:7) = strain_12';
    stress(i,1) = index;
    stress(i,2:4) = stress_xy';
    stress(i,5:7) = stress_12';
end
%% Part 4: Calculate Layer Stresses
% Input:
%   1. Table with layer Sbar,Qbar,T
%   2. Table with layer strains (x,y,xy)
% Output:
%   1. Table with layer stresses (x,y,xy)
%% Part 5: Calculate Layer Principal Stresses, Strains
% Input:
%   1. Table with layer Sbar,Qbar,T
%   2. Table with layer strains (x,y,xy)
%   3. Table with layer stresses (x,y,xy)
% Output:
%   1. Table with layer principal strains
%   2. Table with layer principal stresses
