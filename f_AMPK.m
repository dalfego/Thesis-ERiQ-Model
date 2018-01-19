function AMPK = f_AMPK(ATPr)

global AMPK_SA

AMPK = AMPK_SA*(1./ATPr);

end