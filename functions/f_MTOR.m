function [MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay)

global MTOR_SA

MTORs = (AKT-(AMPK*4));
MTORa = Ay-(MTORs*1.5);

MTOR = MTOR_SA*(1+MTORs+MTORa);

end
