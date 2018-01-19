function GLYCOL = f_GLYCOL(u3,SIRT)

global GLYCOL_SA

gain3 = 0.25;

GLYCOL = GLYCOL_SA*((gain3*u3)+(1./SIRT));

end