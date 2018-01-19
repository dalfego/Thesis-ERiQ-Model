xx = 600;

%% ODE Solving

t0 = [0,xx];

% Conditions to establish homeostais: Y0 = [0 3 1 2 0 0 0 0 1 0 0];
%   Order: [MDAMAGE MFUNCT MENZY GLYCOL GLYENZ Cy Ay Cx Ax Cz Az]

% Homeostatis Initial Conditions:
Y0 = [0.0724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   ...
-0.1936   -0.0000    0.8734   -0.7944    0.0794];
%   For pulse simulation, these initial conditions can be changed ...
%   at specified time to indicate sudden increase in MDAMAGE or MFUNCT

% End simulation when MFUNCT=0.5
options = odeset('Events',@ERiQ_event);

% Modified Damage and p53 rates using global constants across functions
%   Can change to obtain values for 3d multiple parameter plots
global P53_Base P53_Act MDR
P53_Base = 4;
P53_Act = 1;
MDR = 1.8E-3;

% Global used for constants only - for quick sensitivity analysis
%   i.e. change global constant from 1 to 1.1 for 10% increase in node
global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2

FREERAD_SA = 1; GLYCOL_SA = 1; MDAMAGE_SA = 1;

PTEN_SA = 1; AKT_SA = 1; NFKB_SA = 1; P53_SA = 1; 
AMPK_SA = 1; PGC1a_SA = 1; MTOR_SA = 1; AUTO_SA = 1; FOXO_SA = 1;
SIRT_SA = 1; HIF_SA = 1; GLU_SA = 1; PYR_SA = 1; ROS_SA = 1;
MDR_SA = 1; NADr_SA = 1; MENZY_SA = 1; GLYENZ_SA = 1;

ROS_SA2 = 1;

[t,Y] = ode15s(@ERiQ_NFKB_aging,t0,Y0,options);

tnew = [xx 1000000];
Ynew = [1.9535    1.8710   -2.1378    2.9068   -3.0764    0.0083    0.0231   -0.0038    0.7868   -1.1184    0.1118];
NFKB_SA = 0.9;
[t2,Y2] = ode15s(@ERiQ_NFKB_aging,tnew,Ynew,options);

MDAMAGE = Y(:,1); MDAMAGE2 = Y2(:,1); 
MFUNCT = Y(:,2); MFUNCT2 = Y2(:,2);
MENZY = Y(:,3); MENZY2 = Y2(:,3); 
GLYCOL = Y(:,4); GLYCOL2 = Y2(:,4); 
GLYENZ = Y(:,5); GLYENZ2 = Y2(:,5);
Cy = Y(:,6); Ay = Y(:,7); Cy2 = Y2(:,6); Ay2 = Y2(:,7);
Cx = Y(:,8); Ax = Y(:,9); Cx2 = Y2(:,8); Ax2 = Y2(:,9);
Cz = Y(:,10); Az = Y(:,11); Cz2 = Y2(:,10); Az2 = Y2(:,11);

%% Nodes for analysis
 
[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);
[ATPm2,ATPg2,ATPr2] = f_ATP(GLYCOL2, MFUNCT2);
ROS = f_ROS(Az); 
ROS2 = f_ROS(Az2); 
PTEN = f_PTEN(MFUNCT); 
PTEN2 = f_PTEN(MFUNCT2); 
AKT = f_AKT(PTEN,ROS);
AKT2 = f_AKT(PTEN2,ROS2);
AMPK = f_AMPK(ATPr); 
AMPK2 = f_AMPK(ATPr2);
PYR = f_PYR(GLYCOL);
PYR2 = f_PYR(GLYCOL2);
FOXO = f_FOXO(AKT);
FOXO2 = f_FOXO(AKT2);
[MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay);
[MTORs2,MTORa2,MTOR2] = f_MTOR(AKT2,AMPK2,Ay2);
NFKB = f_NFKB_aging(AKT,ROS,MTOR,t);
NFKB2 = f_NFKB_aging(AKT2,ROS2,MTOR2,t2);
NADr = f_NADr(MFUNCT);
NADr2 = f_NADr(MFUNCT2);
[P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax);
[P53s2,P53a2,P532] = f_P53(AKT2,NFKB2,ROS2,Ax2);
AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53);
AUTOPHAGY2 = f_AUTO(MTOR2,FOXO2,ROS2,P532);
MD = f_MD(MFUNCT,ROS); 
MD2 = f_MD(MFUNCT2,ROS2); 
GLU = f_GLU(NFKB);
GLU2 = f_GLU(NFKB2);
HIF = f_HIF(AKT);
HIF2 = f_HIF(AKT2);
SIRT = f_SIRT(NADr);
SIRT2 = f_SIRT(NADr2);
PGC1alpha = f_PGC1a(AMPK,SIRT);
PGC1alpha2 = f_PGC1a(AMPK2,SIRT2);
Uz = f_FREERAD(P53,MDAMAGE,Cz,FOXO);
Uz2 = f_FREERAD(P532,MDAMAGE2,Cz2,FOXO2);

duration = t(end);
duration2 = t2(end)

%% Plot Major Nodes

figure(1);
hold on
yyaxis left
plot(t2,ATPm2,'-xb',t2,ATPg2,'-dg','LineWidth',1.5,'MarkerSize',8);
plot(t2,ROS2,'->k','LineWidth',1.5,'MarkerSize',8);
plot(t2,AMPK2,'-oc',t2,NFKB2,'-sm',t2,P532,'-*r','LineWidth',1.5,'MarkerSize',8);
plot(t2,MDAMAGE2,'-<','Color',1/255*[200,200,0],'LineWidth',1.5,'MarkerSize',8);
plot(t2,MTOR2,'-+','Color',[1 .5 0],'LineWidth',1.5,'MarkerSize',8);
plot(t2,AKT2,'-p','Color',1/255*[0,150,87],'LineWidth',1.5,'MarkerSize',8);
xlabel('Model Time','FontSize',14);
ylabel('Quantity','FontSize',14);
yyaxis right
plot(t2,AUTOPHAGY2,'--','Color',1/255*[120,120,120],'LineWidth',1.5);
ylabel('Quantity (Autophagy Only)','FontSize',14);
legend('ATP: Mito.','ATP: Gly.','ROS','AMPK','NF-kB','P53', ...
    'Mito. Damage','mTOR','AKT','Autophagy')
hold off
