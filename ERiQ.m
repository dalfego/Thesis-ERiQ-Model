function dY = ERiQ(t,Y)

%% Rename variables

MDAMAGE = Y(1);
MFUNCT = Y(2);
MENZY = Y(3);
GLYCOL = Y(4);
GLYENZ = Y(5);
Cy = Y(6); Ay = Y(7);
Cx = Y(8); Ax = Y(9);
Cz = Y(10); Az = Y(11);

[ATPm,ATPg,ATPr] = f_ATP(GLYCOL, MFUNCT);

%% ROS Module

ROS = f_ROS(Az);             
% {* x, at > 0.5 runtime > 900, flips at > 5.0 to reduced runtimes}
PTEN = f_PTEN(MFUNCT);
AKT = f_AKT(PTEN,ROS);

%% mTOR

AMPK = f_AMPK(ATPr);
NADr = f_NADr(MFUNCT);       
% { NAD+/NADH ratio}
SIRT = f_SIRT(NADr);
PGC1alpha = f_PGC1a(AMPK,SIRT);
[MTORs,MTORa,MTOR] = f_MTOR(AKT,AMPK,Ay);
NFKB = f_NFKB(AKT,ROS,MTOR); 
% { negative feedback }

%% p53 Autoregulation

%P53Act = 1 in script {< 1 = no activation, modify this between 0.1 and 4};                   
[P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax);

%% ROS cont.

FOXO = f_FOXO(AKT);           
% { FOXO detoxifies, is inhibited by AKT }
Uz = f_FREERAD(P53,MDAMAGE,Cz,FOXO);  
% Free Radicals

%% Autophagy Module

AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53);       
% { Autophagy depends on mTOR and P53 }

%% Mito / Apoptosis

HIF = f_HIF(AKT);
PYR = f_PYR(GLYCOL);                             
% { a function of Glycolysis }
r2 = (PGC1alpha+PYR+P53)-(HIF*0.2)-(NFKB*0.2);     
% { positive feedback }
MD = f_MD(MFUNCT,ROS);                           
% { this rate will be modified slightly up and down }                    
% { <<<<<< constant scales damage accumulation in next equation, range }
d2 = MDAMAGE;
k4 = (r2-d2);
gain2 = 0.05;                     
% { >>> gain determines speed of response, mitochondria respond 
%slower than glycolysis }
u2 = r2+MENZY;

%% Glycolysis

GLU = f_GLU(NFKB);              
% { GLU is Glucose uptake, inhibited by NFkB, should remain positive }
r3 = GLU+(HIF*.01)+(NADr*.01);
u3 = r3+GLYENZ;
k6 = r3;

%% ODEs
global ROS_SA2
ROS = ROS*ROS_SA2;
dY = zeros(11,1);
dY(1) = f_MDAMAGE(MD,AUTOPHAGY);                % Mitochondrial Damage
dY(2) = f_MFUNCT(gain2,u2,SIRT);                % Mitochondrial Function
dY(3) = f_MENZY(ATPm,k4,MENZY);                 % Mitochondrial Enzymes
dY(4) = f_GLYCOL(u3,SIRT);                      % Glycolysis
dY(5) = f_GLYENZ(ATPg,k6,GLYENZ);               % Glycolytic Enzymes
[dY(6),dY(7)] = f_MTOR_feedback(MTORa,Cy);      % mTOR Feedback loop
[dY(8),dY(9)] = f_P53_feedback(P53a,Cx);        % p53 Feedback loop
[dY(10),dY(11)] = f_ROS_feedback(ROS,Cz,Uz);    % ROS Feedback loop

end