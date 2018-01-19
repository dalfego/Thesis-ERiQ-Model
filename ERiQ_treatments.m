t0 = [0,10000];
Y0 = [0.0724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   ...
-0.1936   -0.0000    0.8734   -0.7944    0.0794];

options = odeset('Events',@ERiQ_event);

global P53_Base P53_Act MDR
P53_Base = 4;
P53_Act = 1;
MDR = 1.8E-3;

%% Normal Conditions
global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2
    
ROS_SA = 1; PTEN_SA = 1; FREERAD_SA = 1; AKT_SA = 1; NFKB_SA = 1;
P53_SA = 1; AMPK_SA = 1; PGC1a_SA = 1; MTOR_SA = 1; MDAMAGE_SA = 1;
MDR_SA = 1; PYR_SA = 1; AUTO_SA = 1; FOXO_SA = 1; NADr_SA = 1;
GLYCOL_SA = 1; GLU_SA = 1; HIF_SA = 1; SIRT_SA = 1; MENZY_SA = 1;
GLYENZ_SA = 1; ROS_SA2 = 1;
[t,Y] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE = Y(:,1); 
MFUNCT = Y(:,2); 
GLYCOL = Y(:,4); 

[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);

%% mTOR Inhibition
% Low inhibition concentration = 0.75, high = 0.55.
% Change these to enact different concentrations

MTOR_SA = 0.75;
[t_mtor_inh1,Y_mtor_inh1] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_mtor_inh1 = Y_mtor_inh1(:,1); 
MFUNCT_mtor_inh1 = Y_mtor_inh1(:,2); 
GLYCOL_mtor_inh1 = Y_mtor_inh1(:,4); 
 
[ATPm_mtor_inh1,ATPg_mtor_inh1,ATPr_mtor_inh1] = f_ATP(GLYCOL_mtor_inh1, MFUNCT_mtor_inh1);

MTOR_SA = 0.55;
[t_mtor_inh2,Y_mtor_inh2] = ode15s(@ERiQ,t0,Y0,options);
p_mtor1 = round(percentdiff(t(end),t_mtor_inh1(end)),2);

MDAMAGE_mtor_inh2 = Y_mtor_inh2(:,1); 
MFUNCT_mtor_inh2 = Y_mtor_inh2(:,2); 
GLYCOL_mtor_inh2 = Y_mtor_inh2(:,4); 
 
[ATPm_mtor_inh2,ATPg_mtor_inh2,ATPr_mtor_inh2] = f_ATP(GLYCOL_mtor_inh2, MFUNCT_mtor_inh2);
p_mtor2 = round(percentdiff(t(end),t_mtor_inh2(end)),2);

%% AKT inhibition
% Low inhibition concentration = 0.75, high = 0.55.
% Change these to enact different concentrations

MTOR_SA = 1;
AKT_SA = 0.75;
[t_akt_inh1,Y_akt_inh1] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_akt_inh1 = Y_akt_inh1(:,1); 
MFUNCT_akt_inh1 = Y_akt_inh1(:,2); 
GLYCOL_akt_inh1 = Y_akt_inh1(:,4); 
 
[ATPm_akt_inh1,ATPg_akt_inh1,ATPr_akt_inh1] = f_ATP(GLYCOL_akt_inh1, MFUNCT_akt_inh1);
p_akt1 = round(percentdiff(t(end),t_akt_inh1(end)),2);

MTOR_SA = 1;
AKT_SA = 0.55;
[t_akt_inh2,Y_akt_inh2] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_akt_inh2 = Y_akt_inh2(:,1); 
MFUNCT_akt_inh2 = Y_akt_inh2(:,2); 
GLYCOL_akt_inh2 = Y_akt_inh2(:,4); 
 
[ATPm_akt_inh2,ATPg_akt_inh2,ATPr_akt_inh2] = f_ATP(GLYCOL_akt_inh2, MFUNCT_akt_inh2);
p_akt2 = round(percentdiff(t(end),t_akt_inh2(end)),2);

%% Autophagy inhibition
% Low inhibition concentration = 0.75, high = 0.55.
% Change these to enact different concentrations

MTOR_SA = 1; AKT_SA = 1;
AUTO_SA = 0.75;
[t_auto_inh1,Y_auto_inh1] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_auto_inh1 = Y_auto_inh1(:,1); 
MFUNCT_auto_inh1 = Y_auto_inh1(:,2); 
GLYCOL_auto_inh1 = Y_auto_inh1(:,4); 
 
[ATPm_auto_inh1,ATPg_auto_inh1,ATPr_auto_inh1] = f_ATP(GLYCOL_auto_inh1, MFUNCT_auto_inh1);
p_auto1 = round(percentdiff(t(end),t_auto_inh1(end)),2);

MTOR_SA = 1; AKT_SA = 1;
AUTO_SA = 1.25;
[t_auto_inh2,Y_auto_inh2] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_auto_inh2 = Y_auto_inh2(:,1); 
MFUNCT_auto_inh2 = Y_auto_inh2(:,2); 
GLYCOL_auto_inh2 = Y_auto_inh2(:,4); 
 
[ATPm_auto_inh2,ATPg_auto_inh2,ATPr_auto_inh2] = f_ATP(GLYCOL_auto_inh2, MFUNCT_auto_inh2);
p_auto2 = round(percentdiff(t(end),t_auto_inh2(end)),2);

%% NFKB inhibition
% Low inhibition concentration = 0.75, high = 0.55.
% Change these to enact different concentrations

MTOR_SA = 1; AKT_SA = 1; AUTO_SA = 1;
NFKB_SA = 0.75;
[t_nfkb_inh1,Y_nfkb_inh1] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_nfkb_inh1 = Y_nfkb_inh1(:,1); 
MFUNCT_nfkb_inh1 = Y_nfkb_inh1(:,2); 
GLYCOL_nfkb_inh1 = Y_nfkb_inh1(:,4); 
 
[ATPm_nfkb_inh1,ATPg_nfkb_inh1,ATPr_nfkb_inh1] = f_ATP(GLYCOL_nfkb_inh1, MFUNCT_nfkb_inh1);
p_nfkb1 = round(percentdiff(t(end),t_nfkb_inh1(end)),2);

MTOR_SA = 1; AKT_SA = 1; AUTO_SA = 1.25;
NFKB_SA = 0.55;
[t_nfkb_inh2,Y_nfkb_inh2] = ode15s(@ERiQ,t0,Y0,options);

MDAMAGE_nfkb_inh2 = Y_nfkb_inh2(:,1); 
MFUNCT_nfkb_inh2 = Y_nfkb_inh2(:,2); 
GLYCOL_nfkb_inh2 = Y_nfkb_inh2(:,4); 
 
[ATPm_nfkb_inh2,ATPg_nfkb_inh2,ATPr_nfkb_inh2] = f_ATP(GLYCOL_nfkb_inh2, MFUNCT_nfkb_inh2);
p_nfkb2 = round(percentdiff(t(end),t_nfkb_inh2(end)),2);

%% 
pp = [ num2str(0),'%'];

figure(1);
hold on
plot(t,ATPm,'b',t_mtor_inh1,ATPm_mtor_inh1,'b--',t_mtor_inh2,ATPm_mtor_inh2,'b-.','LineWidth',2)
plot(t,MDAMAGE,'Color',1/255*[200,200,0],'LineWidth',2);
plot(t_mtor_inh1,MDAMAGE_mtor_inh1,'Color',1/255*[200,200,0],'LineStyle','--','LineWidth',2);
plot(t_mtor_inh2,MDAMAGE_mtor_inh2,'Color',1/255*[200,200,0],'LineStyle','-.','LineWidth',2);
xlabel('Time','FontSize',14);
ylabel('Quantity','FontSize',14);
title('mTOR Inhibition')
legend('ATPm','ATPm, low mTOR inh','ATPm, high mTOR inh','MDamage','MDamage, low mTOR inh','MDamage, high mTOR inh')
p1 = [ num2str(p_mtor1),'%'];
k = vline(t_mtor_inh1(end),'k:',p1);
p2 = [ num2str(p_mtor2),'%'];
k2 = vline(t_mtor_inh2(end),'k:',p2);
k3 = vline(t(end),'k:',pp);
hold off

figure(2);
hold on
plot(t,ATPm,'b',t_akt_inh1,ATPm_akt_inh1,'b--',t_akt_inh2,ATPm_akt_inh2,'b-.','LineWidth',2)
plot(t,MDAMAGE,'Color',1/255*[200,200,0],'LineWidth',2);
plot(t_akt_inh1,MDAMAGE_akt_inh1,'Color',1/255*[200,200,0],'LineStyle','--','LineWidth',2);
plot(t_akt_inh2,MDAMAGE_akt_inh2,'Color',1/255*[200,200,0],'LineStyle','-.','LineWidth',2);
xlabel('Time','FontSize',14);
ylabel('Quantity','FontSize',14);
title('AKT Inhibition')
legend('ATPm','ATPm, low AKT inh','ATPm, high AKT inh','MDamage','MDamage, low AKT inh','MDamage, high AKT inh')
p3 = [ num2str(p_akt1),'%'];
j = vline(t_akt_inh1(end),'k:',p3);
p4 = [ num2str(p_akt2),'%'];
j2 = vline(t_akt_inh2(end),'k:',p4);
j3 = vline(t(end),'k:',pp);
hold off

figure(3);
hold on
plot(t,ATPm,'b',t_auto_inh1,ATPm_auto_inh1,'b--',t_auto_inh2,ATPm_auto_inh2,'b-.','LineWidth',2)
plot(t,MDAMAGE,'Color',1/255*[200,200,0],'LineWidth',2);
plot(t_auto_inh1,MDAMAGE_auto_inh1,'Color',1/255*[200,200,0],'LineStyle','--','LineWidth',2);
plot(t_auto_inh2,MDAMAGE_auto_inh2,'Color',1/255*[200,200,0],'LineStyle','-.','LineWidth',2);
xlabel('Time','FontSize',14);
ylabel('Quantity','FontSize',14);
title('Autophagy Inhibition/Activation')
legend('ATPm','ATPm, Autoph inh','ATPm, Autoph act','MDamage','MDamage, Autoph inh','MDamage, Autoph act')
p5 = [ num2str(p_auto1),'%'];
i = vline(t_auto_inh1(end),'k:',p5);
p6 = [ num2str(p_auto2),'%'];
i2 = vline(t_auto_inh2(end),'k:',p6);
i3 = vline(t(end),'k:',pp);
hold off

figure(4);
hold on
plot(t,ATPm,'b',t_nfkb_inh1,ATPm_nfkb_inh1,'b--',t_nfkb_inh2,ATPm_nfkb_inh2,'b-.','LineWidth',2)
plot(t,MDAMAGE,'Color',1/255*[200,200,0],'LineWidth',2);
plot(t_nfkb_inh1,MDAMAGE_nfkb_inh1,'Color',1/255*[200,200,0],'LineStyle','--','LineWidth',2);
plot(t_nfkb_inh2,MDAMAGE_nfkb_inh2,'Color',1/255*[200,200,0],'LineStyle','-.','LineWidth',2);
xlabel('Time','FontSize',14);
ylabel('Quantity','FontSize',14);
title('NFkB Inhibition')
legend('ATPm','ATPm, low NFkB inh','ATPm, high NFkB inh','MDamage','MDamage, low NFkB inh','MDamage, high NFkB inh')
p7 = [ num2str(p_nfkb1),'%'];
h = vline(t_nfkb_inh1(end),'k:',p7);
p8 = [ num2str(p_auto2),'%'];
h2 = vline(t_nfkb_inh2(end),'k:',p8);
h3 = vline(t(end),'k:',pp);
hold off