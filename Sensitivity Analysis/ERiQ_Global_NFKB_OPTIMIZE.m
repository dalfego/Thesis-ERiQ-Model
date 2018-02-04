
%% ODE Solving

t0 = [250,10000000];

Y0 = [1.0362    2.9400   -1.6770    2.5340   -2.4089    0.0025   -0.1427   -0.0011    0.8624   -0.9848    0.0985];

options = odeset('Events',@ERiQ_event);

global P53_Base P53_Act MDR
P53_Base = 4;
P53_Act = 1;
MDR = 1.8E-3;

global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2
    
FREERAD_SA = 1; MDAMAGE_SA = 1; GLYCOL_SA = 1; ROS_SA = 1;
GLU_SA = 1; PYR_SA = 1; ROS_SA2 = 1; MDR_SA = 1; NADr_SA = 1;

% CHOOSE random values for parameters from Gaussian distribution
% normrnd(1,.3); mean = 1 and std = .3
% Chose 0.03 because mean+3*std = 90% of the data, meaning the values will
% never go above 2 or below 0 (roughly)

PTEN_SA = normrnd(0.95,.0165);   %Decrease in PTEN -> increased life, so random variables decreasing PTEN from 1-10%
AKT_SA = normrnd(0.95,.0165);    %Decrease in AKT -> increased life, so random variables decreasing AKT from 1-10%
NFKB_SA = normrnd(0.95,.0165);   %Decrease in NFKB -> increased life, so random variables decreasing NFKB from 1-10%
P53_SA = normrnd(1.05,.0165);    %Increase in P53 -> increased life, so random variables increasing P53 from 1-10%
AMPK_SA = normrnd(1.05,.0165);   %Increase in AMPK -> increased life, so random variables increasing AMPK from 1-10%
PGC1a_SA = normrnd(0.95,.0165);  %Decrease in PGC1 -> increased life, so random variables decreasing PGC1 from 1-10%
MTOR_SA = normrnd(0.95,.0165);   %Decrease in MTOR -> increased life, so random variables decreasing MTOR from 1-10%
AUTO_SA = normrnd(1.05,.0165);   %Increase in AUTO -> increased life, so random variables increasing AUTO from 1-10%
FOXO_SA = normrnd(1.05,.0165);   %Increase in FOXO -> increased life, so random variables increasing FOXO from 1-10%
SIRT_SA = normrnd(1.05,.0165);   %Increase in SIRT -> increased life, so random variables increasing SIRT from 1-10%
HIF_SA = normrnd(0.95,.0165);    %Decrease in HIF -> increased life, so random variables decreasing HIF from 1-10%
MENZY_SA = normrnd(0.95,.0165);  %Decrease in MENZY -> increased life, so random variables decreasing MENZY from 1-10%
GLYENZ_SA = normrnd(0.95,.0165); %Decrease in GLYENZ -> increased life, so random variables decreasing GLYENZ from 1-10%

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

duration = t(end);

% %% Plot Major Nodes
% 
% figure(1);
% hold on
% plot(t,ATPm,'-xb',t,ATPg,'-dg','LineWidth',1.5,'MarkerSize',8);
% plot(t,ROS,'->k','LineWidth',1.5,'MarkerSize',8);
% plot(t,AMPK,'-oc',t,NFKB,'-sm',t,P53,'-*r','LineWidth',1.5,'MarkerSize',8);
% plot(t,MDAMAGE,'-<','Color',1/255*[200,200,0],'LineWidth',1.5,'MarkerSize',8);
% plot(t,MTOR,'-+','Color',[1 .5 0],'LineWidth',1.5,'MarkerSize',8);
% plot(t,AKT,'-p','Color',1/255*[0,150,87],'LineWidth',1.5,'MarkerSize',8);
% xlabel('Time','FontSize',14);
% ylabel('Quantity','FontSize',14);
% legend('ATP: Mito.','ATP: Gly.','ROS','AMPK','NFKB','P53', ...
%     'Mito. Damage','mTOR','AKT')
% hold off


