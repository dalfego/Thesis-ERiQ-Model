function [ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT)

ATPg = GLYCOL;
ATPm = MFUNCT;
ATPr = ATPg+ATPm;

end