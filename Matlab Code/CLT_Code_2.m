% CLT Code 2.0
% Author(s): Daniel Donate-Perez
% Date: 3/19/2021
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
%% Laminate Input Prompt
% Number of Layers
number_layers = input('Number of Layers: ')

% Layer Material Property Input
matl_prompt = input('Does each layer have the same material properties? (1 = Yes, 0 = No): ')

count = 1;

laminate_properties(:,1) = 1:number_layers;
if matl_prompt == 1
    fprintf('For all layers 1 to %.0f: \n',number_layers);
    laminate_properties(:,2) = input('E1 [Pa]: ');
    laminate_properties(:,3) = input('E2 [Pa]: ');
    laminate_properties(:,4) = input('v12: ');
    laminate_properties(:,5) = input('G12 [Pa]: ');
elseif matl_prompt == 0
    while count <= number_layers
        fprintf('For layer %.0f: \n',count);
        laminate_properties(count,2) = input('E1 [Pa]: ');
        laminate_properties(count,3) = input('E2 [Pa]: ');
        laminate_properties(count,4) = input('v12: ');
        laminate_properties(count,5) = input('G12 [Pa]: ');
        count = count + 1;
    end
else
    return
end

% CTE and CMD input
CTE_prompt = input('Do you want to account for thermal and moisture effects? (1=Yes, 0=No): ')
count = 1;

if CTE_prompt == 1
    if matl_prompt == 1
        fprintf('For all layers 1 to %.0f: \n',number_layers);
        laminate_properties(:,6) = input('CTE 1 [1/Degrees Celsius]: ');
        laminate_properties(:,7) = input('CTE 2 [1/Degrees Celsius]: ');
        laminate_properties(:,8) = input('CTE 3 [1/Degrees Celsius]: ');
        laminate_properties(:,9) = input('CMD 1 [1/% Moisture]: ');
        laminate_properties(:,10) = input('CMD 2 [1/% Moisture]: ');
        laminate_properties(:,11) = input('CMD 3 [1/% Moisture]: ');
    elseif matl_prompt == 0
        while count <= number_layers
            fprintf('For layer %.0f: \n',count);
            laminate_properties(count,6) = input('CTE 1 [1/Degrees Celsius]: ');
            laminate_properties(count,7) = input('CTE 2 [1/Degrees Celsius]: ');
            laminate_properties(count,8) = input('CTE 3 [1/Degrees Celsius]: ');
            laminate_properties(count,9) = input('CMD 1 [% Moisture]: ');
            laminate_properties(count,10) = input('CMD 2 [% Moisture]: ');
            laminate_properties(count,11) = input('CMD 3 [% Moisture]: ');
            count = count + 1;
        end
    else
        return
    end
elseif CTE_prompt == 0
else
    return
end

count = 1;

while count <= number_layers
    fprintf('Orientation of layer %.0f: \n',count);
    laminate_properties(count,12) = input('Layer Orientation [degrees]: ');
    count = count + 1;
end
% Layer Thickness Input
layer_prompt = input('Is each layer the same thickness? (1=Yes, 0=No): ')
count = 1;

if layer_prompt == 1
    fprintf('Thickness of all layers 1 to %.0f: \n',number_layers);
    laminate_properties(:,13) = input('Layer Thickness [m]: ');
    z_vect(1) = -laminate_properties(1,13)*number_layers/2;
    while count <= number_layers
        z_vect(count + 1) = z_vect(1) + laminate_properties(1,13)*count;
        count = count + 1;
    end
elseif layer_prompt == 0
    while count <= number_layers
        fprintf('Thickness of layer %.0f: \n',count);
        laminate_properties(count,13) = input('Layer Thickness [m]: ');
        count = count + 1;
    end
    z_vect(1) = -sum(laminate_properties(:,13))/2
    count = 1;
    while count <= number_layers
        z_vect(count+1) = z_vect(1) + sum(laminate_properties(1:count,13));
        count = count + 1;
    end
end

% Displaying Laminate Matrix
colNames = {'Layer_Number','E1','E2','v12','G12',...
    'CTE_1','CTE_2','CTE_3','CMD_1','CMD_2','CMD_3',...
    'Layer_Orientation','Layer_Thickness'};
laminate_table = array2table(laminate_properties,'VariableNames',colNames)

check_prompt = input('Are these values correct? (1=Yes, 0=No): ')

if check_prompt == 1
else
    return
end

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
%% Part 3: Calculate Layer Principal Stresses and Strains
% Input:
%   1. Force-moment resultants
%   2. abd matrix
%   3. Table with layer S,Qbar,T
% Pre-output:
%   Neutral plane strains and Neutral plane curvatures
% Output
%   1. Table with layer strains (x,y,xy,1,2,3)
%   2. Table with layer stresses (x,y,xy,1,2,3)

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
colNames = {'Layer_Number','strain_x','strain_y',...
            'strain_xy','strain_1','strain_2',...
            'strain_12'};
strain_table = array2table((strain),'VariableNames',colNames)    
colNames = {'Layer_Number','stress_x','stress_y',...
            'stress_xy','stress_1','stress_2',...
            'stress_12'};
stress_table = array2table((stress),'VariableNames',colNames)
stress1T_max = min(stress_table.stress_1(stress_table.stress_1>0))
stress1C_max = max(stress_table.stress_1(stress_table.stress_1<0))
stress2T_max = min(stress_table.stress_2(stress_table.stress_2>0))
stress2C_max = max(stress_table.stress_2(stress_table.stress_2<0))
stress12_max = min(abs(stress_table.stress_12))