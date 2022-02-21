function [NM] = NM(x,y,FM)
NM(1) = FM(1)/y;NM(2) = FM(2)/x;NM(3) = FM(1)/y;
NM(4) = FM(4)/y;NM(5) = -FM(5)/x;NM(6) = -FM(6)/y;
NM = NM'
end

