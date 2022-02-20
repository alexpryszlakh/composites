# CLT code
# Author: Alex Hernandez-Pryszlak, Daniel Donate-Perez
# Date:   2/15/2022
# ---------------------------------------------------------------
"""%% Part 0: Material Properties (Done)
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
number_layers = input('Number of Layers: ');

% Layer Material Property Input
matl_prompt = input('Does each layer have the same material properties? (1 = Yes, 0 = No): '); """
# ----------------------------------------------------------------------------------------------------
import numpy as np

number_layers = int(input('Number of Layers: '))
layers = number_layers + 1

# Layer Material Property Input
mat_prompt = input('Does each layer have the same material properties? (1 = Yes, 0 = No): ')
count = 1
laminate_properties = np.array(range(1, layers))

E1_array = np.array([])
E2_array = np.array([])
v12_array = np.array([])
G12_array = np.array([])

if int(mat_prompt) == 1:
    print('For all layers 1 to ' + str(number_layers) + ':')
    E1 = input('E1 [Pa]: ')
    E1_array = np.full((layers, 1), E1)
    E2 = input('E2 [Pa]: ')
    E2_array = np.full((layers, 1), E2)
    v12 = input('v12: ')
    v12_array = np.full((layers, 1), v12)
    G12 = input('G12 [Pa]: ')
    G12_array = np.full((layers, 1), G12)
else:
    while count <= number_layers:
        print('For layer ' + str(count) + ':')
        E1 = input('E1 [Pa]: ')
        E1_array = np.append(E1_array, E1)
        E2 = input('E2 [Pa]: ')
        E2_array = np.append(E2_array, E2)
        v12 = input('v12: ')
        v12_array = np.append(v12_array, v12)
        G12 = input('G12 [Pa]: ')
        G12_array = np.append(G12_array, G12)
        count = count + 1

count = 1

Lo_array = np.array([])
while count <= number_layers:
    print('Orientation of layer ' + str(count) + ':')
    Lo = input('Layer Orientation [degrees]: ')
    Lo_array = np.append(Lo_array, Lo)
    count = count + 1

# Layer Thickness Input
# --------------------------------------------------------------------------------------------------------
layer_prompt = input('Is each layer the same thickness? (1=Yes, 0=No): ')
count = 1
thickness_array = np.array([])

if layer_prompt == 1:
    print('Thickness of all layers 1 to ' + number_layers + ':')
    thickness = input('Layer Thickness [m]: ')
    thickness_array = np.append(thickness_array, thickness)
    z_vector = -thickness_array * number_layers / 2
    while count <= number_layers:
        z_vect(count + 1) = z_vect(1) + laminate_properties(1, 7) * count
        count = count + 1
elif layer_prompt == 0:
    while count <= number_layers:
        print('Thickness of layer ' + count + ':')
        laminate_properties(count, 7) = input('Layer Thickness [m]: ')
        count = count + 1
    z_vect(1) = -sum(laminate_properties(:, 7)) / 2
    count = 1
    while count <= number_layers:
        z_vect(count + 1) = z_vect(1) + sum(laminate_properties(1:count, 7))
        count = count + 1
# ----------------------------------------------------------------------------------------------------------
'''Make into Pandas DF'''
# Displaying Laminate Matrix
colNames = {'Layer_Number', 'E1', 'E2', 'v12', 'G12', ...,
            'Layer_Orientation', 'Layer_Thickness'};
laminate_table = array2table(laminate_properties, 'VariableNames', colNames)

symmetric = input('Is this laminate symmetric? (1=Yes, 0=No): ')
balanced = input('Is this laminate balanced? (1=Yes, 0=No): ')
crossply = input('Is this laminate crossply? (1=Yes, 0=No): ')

check_prompt = input('Are these values correct? (1=Yes, 0=No): ')

if check_prompt == 1:
else:
    return

'''%% Part 1: Calculate Force-Moment Resultants (Done)
% Input:
%   1. Applied loads
%       a. P(x,y),V(xy),M(x,y,xy)
%   2. Laminate geometry
%       a. Laminate width (along y-axis)
%       b. Laminate length (along x-axis)
% Output:
%   1. Force-Moment Resultants'''

# Input.1: Applied loads

FM(1) = input('Applied normal force in x-direction: ');
FM(2) = input('Applied normal force in y-direction: ');
FM(3) = input('Applied shear force in xy-plane: ');
FM(4) = input('Applied moment about y-axis: ');
FM(5) = input('Applied moment about x-axis: ');
FM(6) = input('Applied twisting moment about xy-plane: ');
# Input.2: Laminate geometry
y = input('Plate width (along y-axis): ');
x = input('Plate length (along x-axis): ');
# Output: Force-Moment Resultants
NM(1) = FM(1) / y;
NM(2) = FM(2) / x;
NM(3) = FM(3) / y;
NM(4) = FM(4) / y;
NM(5) = -FM(5) / x;
NM(6) = -FM(6) / y;
NM = NM
';

''' Part 2: Calculate ABD Matrix (Done)
% Input:
#  1. Table of layer information
#   2. z-vector
# Pre-output:
#   1. Table with Qbar
# Output:
#   1. ABD Matrix
#   2. abd Matrix '''

# Pre-output: Qbar
Qb(1: number_layers, 1) = 1: number_layers;
for i = 1:number_layers:
E1 = laminate_table.E1(i);
E2 = laminate_table.E2(i);
v12 = laminate_table.v12(i);
G12 = laminate_table.G12(i);
theta = laminate_table.Layer_Orientation(i);
Qb(i, 2: 7) = Qbar(E1, E2, v12, G12, theta);

# Output.1: ABD Matrix
ABDmat = ABD(Qb(:, 2: 7), z_vect);
# If symmetric, Bij = 0
if symmetric == 1
    ABDmat(1: 3, 4: 6) = 0;
    ABDmat(4: 6, 1: 3) = 0;

# If balanced, A16 = A26 = 0
if balanced == 1:
    ABDmat(1: 2, 3) = 0;
    ABDmat(3, 1: 2) = 0

# If cross-ply, A16 = A26 = B16 = B26 = D16 = D26 = 0
if crossply == 1
    ABDmat(1: 2, 3) = 0; % A16, A26
    ABDmat(3, 1: 2) = 0; %
    ABDmat(1: 2, 6) = 0; % B16, B26
    ABDmat(6, 1: 2) = 0; %
    ABDmat(4: 5, 3) = 0; %
    ABDmat(3, 4: 5) = 0; %
    ABDmat(4: 5, 6) = 0; % D16, D26
    ABDmat(6, 4: 5) = 0; %
    end

    ''''% Output.2: abd Matrix
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
    %   2. Table with layer stresses (x,y,xy,1,2,3)'''

    # Pre-output:
    e_k = ABDmat\NM
    strain = zeros(2 * number_layers, 7)
    stress = zeros(2 * number_layers, 7)
    zvect_layer = zeros(2 + 2 * (length(z_vect) - 2), 1)
    zvect_layer(1) = z_vect(1)
    zvect_layer(end) = z_vect(end)
    for i = 3:length(zvect_layer) - 1
    index = fix(i / 2 + 0.5)
    zvect_layer(i - 1: i) = z_vect(index)

for i = 1:2 * number_layers
index = fix(i / 2 + 0.5)
strain(i, 1) = index
ex = e_k(1) + zvect_layer(i) * e_k(4)
ey = e_k(2) + zvect_layer(i) * e_k(5)
exy = e_k(3) + zvect_layer(i) * e_k(6)
strain(i, 2: 4) = [ex ey exy]
evect = [ex ey exy]
'
theta = laminate_table.Layer_Orientation(index)
Qb1 = Qb(:, 2: 7)
[stress_xy, stress_12, strain_12] = stress_strain12(Qb1(index,:), evect, theta)
strain(i, 5: 7) = strain_12
'
stress(i, 1) = index
stress(i, 2: 4) = stress_xy
'
stress(i, 5: 7) = stress_12
'

'''Graphing'''
# colNames = {'Layer_Number', 'strain_x', 'strain_y', ...
#             'strain_xy', 'strain_1', 'strain_2', ...
#             'strain_12'}
# strain_table = array2table((strain), 'VariableNames', colNames)
# colNames = {'Layer_Number', 'stress_x', 'stress_y', ...
#             'stress_xy', 'stress_1', 'stress_2', ...
#             'stress_12'};
# stress_table = array2table((stress), 'VariableNames', colNames)
# % Find
# max
# layer
# tensile, compressive, shear
# stress
# stressT1_max = max(stress(:, 4));
# % %

