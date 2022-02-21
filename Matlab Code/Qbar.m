function [Qbar] = Qbar(E1,E2,v12,G12,theta)
% This function calculates the transformed reduced stiffness of a lamina
% The inputs are:
%       E1 = Extensional modulus in longitudinal direction [ Pa ]
%       E2 = Extensional modulus in transverse direction [ Pa ]
%       v12 = Poisson ratio in 1-2 plane
%       G12 = Shear modulus in 1-2 plane [ Pa ]
%       theta = lamina fiber direction [ degrees ], 
%       where:
%           1. +0-degree is longitudinal direction (1-direction)
%           2. +90-degree is transverse direction (2-direction)
%           3. CCW is deemed positive
%       Note: This angle indicates the fiber alignment within the lamina...
%             relative to the coordinate system defined above
%% Calculation of v21 using reciprocity relations ( v21 )
v21 = (E2/E1)*v12; % eqn(2.42)
%% Calculation of reduced stiffness matrix ( [Q] )
Q11 = E1/(1 - v12*v21); % eqn(4.17), #1
Q12 = (v12*E2)/(1 - v12*v21); % eqn(4.17), #2
Q22 = E2/(1 - v12*v21); % eqn(4.17), #3
Q66 = G12; % eqn(4.17), #4
%% Calculation of transformed reduced stiffness matrix
m = cosd(theta);
n = sind(theta);
Qb11 = Q11*m^4 + 2*(Q12 + 2*Q66)*n^2*m^2 + Q22*n^4; % eqn(5.84), #1
Qb12 = (Q11 + Q22 - 4*Q66)*n^2*m^2 + Q12*(n^4 + m^4); % eqn(5.84), #2
Qb16 = (Q11 - Q12 - 2*Q66)*n*m^3 + (Q12 - Q22 + 2*Q66)*n^3*m; % eqn(5.84), #3
Qb22 = Q11*n^4 + 2*(Q12 + 2*Q66)*n^2*m^2 + Q22*m^4; % eqn(5.84), #4
Qb26 = (Q11 - Q12 - 2*Q66)*n^3*m + (Q12 - Q22 + 2*Q66)*n*m^3; % eqn(5.84), #5
Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*n^2*m^2 + Q66*(n^4 + m^4); % eqn(5.84), #6
Qbar = [Qb11 Qb12 Qb16 Qb22 Qb26 Qb66];
end

