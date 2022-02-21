function [ABD] = ABD(Qbar,z_vect)


K = length(z_vect) - 1; %count number of layers in laminate
z = z_vect;
A = zeros(3);B = zeros(3);D = zeros(3);
for i = 1:K
    n = i + 1;
    A11 = Qbar(i,1)*(z(n) - z(i));
    A12 = Qbar(i,2)*(z(n) - z(i));
    A16 = Qbar(i,3)*(z(n) - z(i));
    A22 = Qbar(i,4)*(z(n) - z(i));
    A26 = Qbar(i,5)*(z(n) - z(i));
    A66 = Qbar(i,6)*(z(n) - z(i));

    B11 = 0.5*Qbar(i,1)*(z(n)^2 - z(i)^2);
    B12 = 0.5*Qbar(i,2)*(z(n)^2 - z(i)^2);
    B16 = 0.5*Qbar(i,3)*(z(n)^2 - z(i)^2);
    B22 = 0.5*Qbar(i,4)*(z(n)^2 - z(i)^2);
    B26 = 0.5*Qbar(i,5)*(z(n)^2 - z(i)^2);
    B66 = 0.5*Qbar(i,6)*(z(n)^2 - z(i)^2);

    D11 = (1/3)*Qbar(i,1)*(z(n)^3 - z(i)^3);
    D12 = (1/3)*Qbar(i,2)*(z(n)^3 - z(i)^3);
    D16 = (1/3)*Qbar(i,3)*(z(n)^3 - z(i)^3);
    D22 = (1/3)*Qbar(i,4)*(z(n)^3 - z(i)^3);
    D26 = (1/3)*Qbar(i,5)*(z(n)^3 - z(i)^3);
    D66 = (1/3)*Qbar(i,6)*(z(n)^3 - z(i)^3);
    
    A = A + [A11,A12,A16;
            A12,A22,A26;
            A16,A26,A66];
    B = B + [B11,B12,B16;
            B12,B22,B26;
            B16,B26,B66];
    D = D + [D11,D12,D16;
            D12,D22,D26;
            D16,D26,D66];
end
ABD = [A,B;
       B,D];     

