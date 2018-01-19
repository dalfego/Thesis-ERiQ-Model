%% ODE Solving

t0 = [0,10000];

Y0 = [0 3 1 2 0 0 0 0 1 0 0];
%Y0 = [0.0724    3.6239   -1.3358    2.4010   -2.1968   -0.0000   ...
%-0.1936   -0.0000    0.8734   -0.7944    0.0794];

options = odeset('Events',@ERiQ_event);

global P53_Base P53_Act MDR
P53_Base = 4;
P53_Act = 1;
MDR = 1.8E-3;

global PTEN_SA AKT_SA FREERAD_SA NFKB_SA P53_SA AMPK_SA PGC1a_SA ... 
    MTOR_SA MDAMAGE_SA MDR_SA PYR_SA AUTO_SA FOXO_SA NADr_SA ... 
    GLYCOL_SA GLU_SA HIF_SA SIRT_SA ROS_SA MENZY_SA GLYENZ_SA ROS_SA2
    
ROS_SA = 1; PTEN_SA = 1; FREERAD_SA = 1; AKT_SA = 1; NFKB_SA = 1;
P53_SA = 1; AMPK_SA = 1; PGC1a_SA = 1; MTOR_SA = 1; MDAMAGE_SA = 1;
MDR_SA = 1; PYR_SA = 1; AUTO_SA = 1; FOXO_SA = 1; NADr_SA = 1;
GLYCOL_SA = 1; GLU_SA = 1; HIF_SA = 1; SIRT_SA = 1; MENZY_SA = 1;
GLYENZ_SA = 1; ROS_SA2 = 1;

[t1,Y1] = ode15s(@ERiQ,t0,Y0,options);
MDR = 2E-3;
[t2,Y2] = ode15s(@ERiQ,t0,Y0,options);
MDR = 2.2E-3;
[t3,Y3] = ode15s(@ERiQ,t0,Y0,options);
MDR = 1.6E-3;
[t4,Y4] = ode15s(@ERiQ,t0,Y0,options);
MDR = 1.4E-3;
[t5,Y5] = ode15s(@ERiQ,t0,Y0,options);
MDR = 0.85E-3;
[t6,Y6] = ode15s(@ERiQ,t0,Y0,options);
MDR = 1E-3;
[t7,Y7] = ode15s(@ERiQ,t0,Y0,options);

MFUNCT1 = Y1(:,2); GLYCOL1 = Y1(:,3);
MFUNCT2 = Y2(:,2); GLYCOL2 = Y2(:,3);
MFUNCT3 = Y3(:,2); GLYCOL3 = Y3(:,3);
MFUNCT4 = Y4(:,2); GLYCOL4 = Y4(:,3);
MFUNCT5 = Y5(:,2); GLYCOL5 = Y5(:,3);
MFUNCT6 = Y6(:,2); GLYCOL6 = Y6(:,3);
MFUNCT7 = Y7(:,2); GLYCOL7 = Y7(:,3);


%% Nodes for analysis
 
[ATPm1,ATPg1,ATPr1] = f_ATP(GLYCOL1, MFUNCT1);
[ATPm2,ATPg2,ATPr2] = f_ATP(GLYCOL2, MFUNCT2);
[ATPm3,ATPg3,ATPr3] = f_ATP(GLYCOL3, MFUNCT3);
[ATPm4,ATPg4,ATPr4] = f_ATP(GLYCOL4, MFUNCT4);
[ATPm5,ATPg5,ATPr5] = f_ATP(GLYCOL5, MFUNCT5);
[ATPm6,ATPg6,ATPr6] = f_ATP(GLYCOL6, MFUNCT6);
[ATPm7,ATPg7,ATPr7] = f_ATP(GLYCOL7, MFUNCT7);

%% Plot Major Nodes

figure(1);

hold on
plot(t3,ATPm3,'c:','LineWidth',2);
plot(t2,ATPm2,'c--','LineWidth',2);
plot(t1,ATPm1,'b','LineWidth',2);
plot(t4,ATPm4,'g--','LineWidth',2);
plot(t5,ATPm5,'g:','LineWidth',2);
plot(t7,ATPm7,'r--','LineWidth',2);
plot(t6,ATPm6,'r:','LineWidth',2);
xlim([0 3100]);
legend('MDR = 2.2E-3','MDR = 2.0E-3','MDR = 1.8E-3','MDR = 1.6E-3','MDR = 1.4E-3','MDR = 1.0E-3','MDR = 0.85E-3');
xlabel('Time','FontSize',14);
ylabel('Quantity','FontSize',14);
title('Mitochondrial Function per Damage Rate');
hold off