%% Homeostatic Simulation

% End simulation when MFUNCT=0.5
options = odeset('Events',@ERiQ_event);

% Modified Damage and p53 rates using global constants across functions
global P53_Base P53_Act MDR;
P53_Base = 4;
P53_Act = 1;
MDR = .85E-3;
%Homestasis: P53Act=1, MDR<1E-3
    
% Global used for constants only - for quick sensitivity analysis
global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2
    
ROS_SA = 1; PTEN_SA = 1; FREERAD_SA = 1; AKT_SA = 1; NFKB_SA = 1;
P53_SA = 1; AMPK_SA = 1; PGC1a_SA = 1; MTOR_SA = 1; MDAMAGE_SA = 1;
MDR_SA = 1; PYR_SA = 1; AUTO_SA = 1; FOXO_SA = 1; NADr_SA = 1;
GLYCOL_SA = 1; GLU_SA = 1; HIF_SA = 1; SIRT_SA = 1; MENZY_SA = 1;
GLYENZ_SA = 1; ROS_SA2 = 1;

% Time to reach homeostasis
t0 = [0,3000];
% Initial conditions
Y0 = [0 3 1 2 0 0 0 0 1 0 0];   
[T,Y] = ode15s(@ERiQ,t0,Y0,options);

% Hometostatic conditions achieved at t=3000: [0.0724    3.6239  ...
% -1.3358    2.4010   -2.1968   -0.0000   -0.1936   -0.0000    ... 
%0.8734   -0.7944    0.0794];

%% Mito Damage Pulse

% Pulse at t=3000
t2 = [3000,6000];
% Increased MDAMAGE by 0.5
y2 = [0.5724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   ... 
-0.1936   -0.0000    0.8734   -0.7944    0.0794];
[T2,Y2] = ode15s(@ERiQ,t2,y2,options);

%% Mito Function Negative Pulse

% Pulse at t=3000
t3 = [3000,6000];
% Decreased MFUNCT by 1
y3 = [0.0724    2.6239   -1.3358    2.4010   -2.1968   -0.0000   ...
-0.1936   -0.0000    0.8734   -0.7944    0.0794];
[T3,Y3] = ode15s(@ERiQ,t3,y3,options);

%% Establishing Nodes @ different conditions
% (1)Homeostasis (2)Mito. Damage Pulse (3)Neg. Mito. Function Pulse

MDAMAGE = Y(:,1);  MDAMAGE2 = Y2(:,1); MDAMAGE3 = Y3(:,1);
MFUNCT = Y(:,2);  MFUNCT2 = Y2(:,2); MFUNCT3 = Y3(:,2);
MENZY = Y(:,3);  MENZY2 = Y2(:,3); MENZY3 = Y3(:,3);
GLYCOL = Y(:,4);  GLYCOL2 = Y2(:,4); GLYCOL3 = Y3(:,4);
GLYENZ = Y(:,5);  GLYENZ2 = Y2(:,5); GLYENZ3 = Y3(:,5);
Cy = Y(:,6);  Cy2 = Y2(:,6); Cy3 = Y3(:,6);
Ay = Y(:,7);  Ay2 = Y2(:,7); Ay3 = Y3(:,7);
Cx = Y(:,8);  Cx2 = Y2(:,8); Cx3 = Y3(:,8);
Ax = Y(:,9);  Ax2 = Y2(:,9); Ax3 = Y3(:,9);
Cz = Y(:,10);  Cz2 = Y2(:,10); Cz3 = Y3(:,10);
Az = Y(:,11);  Az2 = Y2(:,11); Az3 = Y3(:,11);

[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);
[ATPm2,ATPg2,ATPr2] = f_ATP(GLYCOL2, MFUNCT2);
[ATPm3,ATPg3,ATPr3] = f_ATP(GLYCOL3, MFUNCT3);

ROS = f_ROS(Az); ROS2 = f_ROS(Az2); ROS3 = f_ROS(Az3); 

PTEN = f_PTEN(MFUNCT); PTEN2 = f_PTEN(MFUNCT2);
PTEN3 = f_PTEN(MFUNCT3); 

AKT = f_AKT(PTEN,ROS); AKT2 = f_AKT(PTEN2,ROS2);
AKT3 = f_AKT(PTEN3,ROS3);

AMPK = f_AMPK(ATPr); AMPK2 = f_AMPK(ATPr2); AMPK3 = f_AMPK(ATPr3);

PYR = f_PYR(GLYCOL); PYR2 = f_PYR(GLYCOL2); PYR3 = f_PYR(GLYCOL3);

FOXO = f_FOXO(AKT); FOXO2 = f_FOXO(AKT2); FOXO3 = f_FOXO(AKT3);

[MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay);
[MTORs2,MTORa2,MTOR2] = f_MTOR(AKT2,AMPK2,Ay2);
[MTORs3,MTORa3,MTOR3] = f_MTOR(AKT3,AMPK3,Ay3);

NFKB = f_NFKB(AKT,ROS,MTOR); NFKB2 = f_NFKB(AKT2,ROS2,MTOR2);
NFKB3 = f_NFKB(AKT3,ROS3,MTOR3);

NADr = f_NADr(MFUNCT); NADr2 = f_NADr(MFUNCT2);
NADr3 = f_NADr(MFUNCT3);

[P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax);
[P53s2,P53a2,P532] = f_P53(AKT2,NFKB2,ROS2,Ax2);
[P53s3,P53a3,P533] = f_P53(AKT3,NFKB3,ROS3,Ax3);

AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53);
AUTOPHAGY2 = f_AUTO(MTOR2,FOXO2,ROS2,P532);
AUTOPHAGY3 = f_AUTO(MTOR3,FOXO3,ROS3,P533);

MD = f_MD(MFUNCT,ROS); MD2 = f_MD(MFUNCT2,ROS2);
MD3 = f_MD(MFUNCT3,ROS3);

GLU = f_GLU(NFKB); GLU2 = f_GLU(NFKB2); GLU3 = f_GLU(NFKB3);

HIF = f_HIF(AKT); HIF2 = f_HIF(AKT2); HIF3 = f_HIF(AKT3);

SIRT = f_SIRT(NADr); SIRT2 = f_SIRT(NADr2); SIRT3 = f_SIRT(NADr3);

PGC1alpha = f_PGC1a(AMPK,SIRT); PGC1alpha2 = f_PGC1a(AMPK2,SIRT2);
PGC1alpha3 = f_PGC1a(AMPK3,SIRT3);

Uz = f_FREERAD(P53,MDAMAGE,Cz,FOXO);
Uz2 = f_FREERAD(P532,MDAMAGE2,Cz2,FOXO2);
Uz3 = f_FREERAD(P533,MDAMAGE3,Cz3,FOXO3);

duration = T(end); duration2 = T2(end); duration3 = T3(end);

%% plot of Mitochondrial Function Pulse

figure(1)
yyaxis left
hold on
plot(T,ATPm,'-xb','MarkerSize',8,'LineWidth',1.5)
plot(T,ATPg,'-dg','MarkerSize',8,'LineWidth',1.5)
plot(T,NFKB,'-sm','MarkerSize',8,'LineWidth',1.5)
plot(T,P53,'-*r','MarkerSize',8,'LineWidth',1.5)
xlabel('Model Time','FontSize',14);
ylabel('Quantity (ATP: Mito, ATP: Gly, NF-kB, p53)','FontSize',14);
h = vline(3000,'k:','Negative Mitochondrial Function Pulse');
xlim([2000 4000])
yyaxis right
plot(T,AMPK,'-oc','MarkerSize',8,'LineWidth',1.5)
plot(T,AKT,'-p','Color',1/255*[0,150,87],'MarkerSize',8,'LineWidth',1.5)
plot(T,ROS,'->k','MarkerSize',8,'LineWidth',1.5)
plot(T,MTOR,'-+','Color',[1 .5 0],'MarkerSize',8,'LineWidth',1.5)
ylabel('Quantity (AMPK, AKT, ROS, mTOR)','FontSize',14);
title('Negative Mitochondrial Function Pulse Recovery')
legend('ATP: Mito.','ATP: Gly','NF-kB','p53','AMPK','AKT','ROS','mTOR')
xlim([2000 4000])

yyaxis left
plot(T3,ATPm3,'-xb','MarkerSize',3,'LineWidth',1.5)
plot(T3,ATPg3,'-dg','MarkerSize',3,'LineWidth',1.5)
plot(T3,NFKB3,'-sm','MarkerSize',3,'LineWidth',1.5)
plot(T3,P533,'-*r','MarkerSize',3,'LineWidth',1.5)
yyaxis right
plot(T3,AMPK3,'-oc','MarkerSize',3,'LineWidth',1.5)
plot(T3,AKT3,'-p','Color',1/255*[0,150,87],'MarkerSize',3,'LineWidth',1.5)
plot(T3,ROS3,'->k','MarkerSize',3,'LineWidth',1.5)
plot(T3,MTOR3,'-+','Color',[1 .5 0],'MarkerSize',3,'LineWidth',1.5)
hold off

%% plot of Mitochondrial Damage Pulses

figure(2)
yyaxis left
hold on
plot(T,ATPm,'-xb','MarkerSize',8,'LineWidth',1.5)
plot(T,ATPg,'-dg','MarkerSize',8,'LineWidth',1.5)
plot(T,NFKB,'-sm','MarkerSize',8,'LineWidth',1.5)
plot(T,P53,'-*r','MarkerSize',8,'LineWidth',1.5)
xlabel('Model Time','FontSize',14);
ylabel('Quantity (ATP: Mito, ATP: Gly, NF-kB, p53)','FontSize',14);
h = vline(3000,'k:','Mitochondrial Damage Pulse');
xlim([2000 4000])
yyaxis right
plot(T,AMPK,'-oc','MarkerSize',8,'LineWidth',1.5)
plot(T,AKT,'-p','Color',1/255*[0,150,87],'MarkerSize',8,'LineWidth',1.5)
plot(T,ROS,'->k','MarkerSize',8,'LineWidth',1.5)
plot(T,MTOR,'-+','Color',[1 .5 0],'MarkerSize',8,'LineWidth',1.5)
ylabel('Quantity (AMPK, AKT, ROS, mTOR)','FontSize',14);
title('Mitochondrial Damage Pulse Recovery')
legend('ATP: Mito.','ATP: Gly','NF-kB','p53','AMPK','AKT','ROS','mTOR')
xlim([2000 4000])

yyaxis left
plot(T2,ATPm2,'-xb','MarkerSize',3,'LineWidth',1.5)
plot(T2,ATPg2,'-dg','MarkerSize',3,'LineWidth',1.5)
plot(T2,NFKB2,'-sm','MarkerSize',3,'LineWidth',1.5)
plot(T2,P532,'-*r','MarkerSize',3,'LineWidth',1.5)
yyaxis right
plot(T2,AMPK2,'-oc','MarkerSize',3,'LineWidth',1.5)
plot(T2,AKT2,'-p','Color',1/255*[0,150,87],'MarkerSize',3,'LineWidth',1.5);
plot(T2,ROS2,'->k','MarkerSize',3,'LineWidth',1.5)
plot(T2,MTOR2,'-+','Color',[1 .5 0],'MarkerSize',3,'LineWidth',1.5)
hold off