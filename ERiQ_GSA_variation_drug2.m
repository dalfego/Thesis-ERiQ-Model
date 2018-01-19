
%% ODE Solving

t0 = [250,1000000];
Y0 = [MDAMAGE(end) MFUNCT(end) MENZY(end) GLYCOL(end) GLYENZ(end) Cy(end) Ay(end) Cx(end) Ax(end) Cz(end) Az(end)];

options = odeset('Events',@ERiQ_event);

global P53_Base P53_Act MDR
P53_Base = 4;
P53_Act = 1;
MDR = 1.8E-3;

global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2

FREERAD_SA = 1; GLYCOL_SA = 1; MDAMAGE_SA = 1;

PTEN_SA = 0.9; AKT_SA = 0.9; NFKB_SA = 0.9; P53_SA = 1.1; 
AMPK_SA = 1.1; PGC1a_SA = 0.9; MTOR_SA = 0.9; AUTO_SA = 1.1; FOXO_SA = 1.1;
SIRT_SA = 1.1; HIF_SA = 0.9; GLU_SA = 1; PYR_SA = 1; ROS_SA = 1;
MDR_SA = 1; NADr_SA = 1; MENZY_SA = 0.9; GLYENZ_SA = 0.9;
ROS_SA2 = 1;

[t2,Y2] = ode15s(@ERiQ_NFKB_aging,t0,Y0,options);

MDAMAGE2 = Y2(:,1); 
MFUNCT2 = Y2(:,2); 
MENZY2 = Y2(:,3); 
GLYCOL2 = Y2(:,4); 
GLYENZ2 = Y2(:,5); 
Cy2 = Y2(:,6); Ay2 = Y2(:,7);
Cx2 = Y2(:,8); Ax2 = Y2(:,9);
Cz2 = Y2(:,10); Az2 = Y2(:,11);

%% Nodes for analysis
 
[ATPm2,ATPg2,ATPr2] = f_ATP(GLYCOL2, MFUNCT2);
ROS2 = f_ROS(Az2); 
PTEN2 = f_PTEN(MFUNCT2); 
AKT2 = f_AKT(PTEN2,ROS2);  
AMPK2 = f_AMPK(ATPr2); 
PYR2 = f_PYR(GLYCOL2);
FOXO2 = f_FOXO(AKT2);
[MTORs2,MTORa2,MTOR2] = f_MTOR(AKT2,AMPK2,Ay2);
NFKB2 = f_NFKB_aging(AKT2,ROS2,MTOR2,t2);
NADr2 = f_NADr(MFUNCT2);
[P53s2,P53a2,P532] = f_P53(AKT2,NFKB2,ROS2,Ax2);
AUTOPHAGY2 = f_AUTO(MTOR2,FOXO2,ROS2,P532);
MD2 = f_MD(MFUNCT2,ROS2); 
GLU2 = f_GLU(NFKB2);
HIF2 = f_HIF(AKT2);
SIRT2 = f_SIRT(NADr2);
PGC1alpha2 = f_PGC1a(AMPK2,SIRT2);
Uz2 = f_FREERAD(P532,MDAMAGE2,Cz2,FOXO2);

duration2 = t2(end);

%% Plot Major Nodes

% figure(1);
% hold on
% yyaxis left
% plot(t,ATPm,'-xb',t,ATPg,'-dg','LineWidth',1.5,'MarkerSize',8);
% plot(t,ROS,'->k','LineWidth',1.5,'MarkerSize',8);
% plot(t,AMPK,'-oc',t,NFKB,'-sm',t,P53,'-*r','LineWidth',1.5,'MarkerSize',8);
% plot(t,MDAMAGE,'-<','Color',1/255*[200,200,0],'LineWidth',1.5,'MarkerSize',8);
% plot(t,MTOR,'-+','Color',[1 .5 0],'LineWidth',1.5,'MarkerSize',8);
% plot(t,AKT,'-p','Color',1/255*[0,150,87],'LineWidth',1.5,'MarkerSize',8);
% xlabel('Model Time','FontSize',14);
% ylabel('Quantity','FontSize',14);
% yyaxis right
% plot(t,AUTOPHAGY,'--','Color',1/255*[120,120,120],'LineWidth',1.5);
% ylabel('Quantity (Autophagy Only)','FontSize',14);
% legend('ATP: Mito.','ATP: Gly.','ROS','AMPK','NF-kB','P53', ...
%     'Mito. Damage','mTOR','AKT','Autophagy')
% hold off