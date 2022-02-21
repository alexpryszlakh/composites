function [laminaEff] = LamEffElastic(prop_vect)
E1 = prop_vect(1);E2 = prop_vect(2);
G12 = prop_vect(3);v12 = prop_vect(4); theta = prop_vect(5);
% Calculate Lamina Effective Elastic Constants wrt Global Coord. x-y-z
s = sind(theta);c = cosd(theta);
Ex = 1/( c^4/E1 + (1/G12 - 2*v12/E1)*s^2*c^2 + s^4/E2 );
Ey = 1/( s^4/E1 + (1/G12 - 2*v12/E1)*s^2*c^2 + c^4/E2 );
Gxy = 1/( (s^4 + c^4)/G12 + 4*(1/E1 + 1/E2 + 2*v12/E1 - 1/(2*G12))*s^2*c^2 );
vxy = Ex*( v12*(s^4 + c^4)/E1 - (1/E1 + 1/E2 - 1/G12)*s^2*c^2 ); 
laminaEff = [Ex,Ey,Gxy,vxy,theta];
% lamina_elastic = array2table(laminaEff,'VariableNames',...
%                          {'Ex','Ey','Gxy','vxy',...
%                           'Rel_Orient'},...
%                           'RowNames',{'123_Local','xyz_Global'});
end

