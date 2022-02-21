%% Tsai-Wu
F11 = 1/(s_L_t*s_L_c);
F1 = 1/s_L_t - 1/s_L_c;
F22 = 1/(s_T_t*s_T_c);
F2 = 1/s_T_t - 1/s_T_c;
F66 = 1/s_LT^2;
F12 = -sqrt(F11*F22)/2;

F11*stress1^2 + F22*stress2^2 + F66*stress6^2 +...
F1*stress1 + F2*stress2 + 2*F12*stress1*stress2 = 1;
%% Effective Engineering Properties
% In terms of abd
E_x = 1/(a11*H);
E_y = 1/(a22*H);
G_xy = 1/(a66*H);
v_xy = -a12/a11;
v_yx = -a12/a22;
% In terms of ABD
E_x = (A11*A22 - A12^2)/(A22*H);
E_y = (A11*A22 - A12^2)/(A11*H);
G_xy = A66/H;
v_xy = A12/A22;
v_yx = A12/A11;

