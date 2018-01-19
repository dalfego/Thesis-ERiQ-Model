t0 = [0,10000000000];
%Y0 = [0 3 1 1 0 0 0 0 1 0 0];
Y0 = [0.0724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   -0.1936   -0.0000    0.8734   -0.7944    0.0794];

global P53_Base
P53_Base = 4;

global PTEN_SA AKT_SA FREERAD_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA GLYENZ_SA MENZY_SA ROS_SA2
    
ROS_SA = 1; PTEN_SA = 1; FREERAD_SA = 1; AKT_SA = 1;
P53_SA = 1; AMPK_SA = 1; PGC1a_SA = 1; MTOR_SA = 1; MDAMAGE_SA = 1;
MDR_SA = 1; PYR_SA = 1; AUTO_SA = 1; FOXO_SA = 1; NADr_SA = 1;
GLYCOL_SA = 1; GLU_SA = 1; HIF_SA = 1; SIRT_SA = 1; GLYENZ_SA = 1;
MENZY_SA = 1; ROS_SA2 = 1;

options = odeset('Events',@ERiQ_event);
[t,Y] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE = Y(:,1); 
MFUNCT = Y(:,2); 
MENZY = Y(:,3); 
GLYCOL = Y(:,4); 
GLYENZ = Y(:,5); 
Cy = Y(:,6); Ay = Y(:,7);
Cx = Y(:,8); Ax = Y(:,9);
Cz = Y(:,10); Az = Y(:,11);

[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);
ROS = f_ROS(Az); 
PTEN = f_PTEN(MFUNCT); 
AKT = f_AKT(PTEN,ROS);  
AMPK = f_AMPK(ATPr); 
PYR = f_PYR(GLYCOL);
FOXO = f_FOXO(AKT);
[MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay);
NFKB = f_NFKB(AKT,ROS,MTOR);
NADr = f_NADr(MFUNCT);
[P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax);
AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53);
MD = f_MD(MFUNCT,ROS); 
GLU = f_GLU(NFKB);
HIF = f_HIF(AKT);
SIRT = f_SIRT(NADr);
PGC1alpha = f_PGC1a(AMPK,SIRT);
Uz = f_FREERAD(P53,MDAMAGE,Cz,FOXO);

duration = t(end);