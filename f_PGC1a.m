function PGC1alpha = f_PGC1a(AMPK,SIRT)

global PGC1a_SA

PGC1alpha = PGC1a_SA*(AMPK+(SIRT*0.1));

end