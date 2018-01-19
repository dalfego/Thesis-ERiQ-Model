function GLYENZ = f_GLYENZ(ATPg,k6,GLYENZ)

global GLYENZ_SA
k5 = -1;

GLYENZ = GLYENZ_SA*((k5*ATPg)-(k6*GLYENZ));

end