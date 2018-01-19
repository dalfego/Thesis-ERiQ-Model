
%% ODE Solving

t0 = [0,250];
% t0 = [250,100000];

% Conditions to establish homeostais: Y0 = [0 3 1 2 0 0 0 0 1 0 0];
%   Order: [MDAMAGE MFUNCT MENZY GLYCOL GLYENZ Cy Ay Cx Ax Cz Az]

% Homeostatis Initial Conditions:
% Y0 = [1.0362    2.9400   -1.6770    2.5340   -2.4089    0.0025   -0.1427   -0.0011    0.8624   -0.9848    0.0985];
Y0 = [0.0724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   -0.1936   -0.0000    0.8734   -0.7944    0.0794];
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

MDAMAGE = Y(:,1); 
MFUNCT = Y(:,2); 
MENZY = Y(:,3); 
GLYCOL = Y(:,4); 
GLYENZ = Y(:,5); 
Cy = Y(:,6); Ay = Y(:,7);
Cx = Y(:,8); Ax = Y(:,9);
Cz = Y(:,10); Az = Y(:,11);

%% Nodes for analysis
 
[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);
ROS = f_ROS(Az); 
PTEN = f_PTEN(MFUNCT); 
AKT = f_AKT(PTEN,ROS);  
AMPK = f_AMPK(ATPr); 
PYR = f_PYR(GLYCOL);
FOXO = f_FOXO(AKT);
[MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay);
NFKB = f_NFKB_aging(AKT,ROS,MTOR,t);
NADr = f_NADr(MFUNCT);
[P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax);
AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53);
MD = f_MD(MFUNCT,ROS); 
GLU = f_GLU(NFKB);
HIF = f_HIF(AKT);
SIRT = f_SIRT(NADr);
PGC1alpha = f_PGC1a(AMPK,SIRT);
Uz = f_FREERAD(P53,MDAMAGE,Cz,FOXO);
ScaleFactor = ROS - .9;

duration = t(end)

%% Plot Major Nodes

figure(1);
hold on
yyaxis left
plot(t,ATPm,'-xb',t,ATPg,'-dg','LineWidth',1.5,'MarkerSize',8);
plot(t,ROS,'->k','LineWidth',1.5,'MarkerSize',8);
plot(t,AMPK,'-oc',t,NFKB,'-sm',t,P53,'-*r','LineWidth',1.5,'MarkerSize',8);
plot(t,MDAMAGE,'-<','Color',1/255*[200,200,0],'LineWidth',1.5,'MarkerSize',8);
plot(t,MTOR,'-+','Color',[1 .5 0],'LineWidth',1.5,'MarkerSize',8);
plot(t,AKT,'-p','Color',1/255*[0,150,87],'LineWidth',1.5,'MarkerSize',8);
xlabel('Model Time','FontSize',14);
ylabel('Quantity','FontSize',14);
yyaxis right
plot(t,AUTOPHAGY,'--','Color',1/255*[120,120,120],'LineWidth',1.5);
ylabel('Quantity (Autophagy Only)','FontSize',14);
legend('ATP: Mito.','ATP: Gly.','ROS','AMPK','NF-kB','P53', ...
    'Mito. Damage','mTOR','AKT','Autophagy')
hold off
%% Parameter Table
%Sensitivity Analysis - node parameter perturbations, can change simulation 
%ending time to evaluate MFUNCT (i.e. at t=500) or run until end to see 
%how perturbation changes lifespan. Print out final values for each node 
%in table

mfunctend = repelem(MFUNCT(end),18);
life = repelem(t(end),18);
X = [PTEN(end) ROS(end) AKT(end) NFKB(end) P53(end) AMPK(end) ...
    PGC1alpha(end) MTOR(end) MDAMAGE(end) MD(end) PYR(end) ...
    AUTOPHAGY(end) FOXO(end) NADr(end) GLYCOL(end) GLU(end) HIF(end) ...
    SIRT(end); mfunctend; life]';
colNames = {'Parameter','MFUNCT','Life'};
rowNames = {'PTEN','FREERAD','AKT','NFKB','P53','AMPK','PGC1a','MTOR', ...
    'MDam','MD','PYR','AUTO','FOXO','NADr','GLYCOL','GLU','HIF','SIRT'};
Parameters = array2table(X,'RowNames',rowNames,'VariableNames',colNames);