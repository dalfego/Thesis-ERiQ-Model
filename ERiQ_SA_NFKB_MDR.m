%% Setup
global P53_Act NFKB_SA MDR ROS_SA2
ROS_SA2 = 1;
P53_Act = 1;

% Set range of desired parameters
NFKB_param = 0.2:0.2:1.8;
damrate = (1.8E-3):(0.1E-3):(2.6E-3);

n = numel(NFKB_param);
m = numel(damrate);

x = zeros(n,1);
nfkb = x;
mdr = x;
life = x;

%% NFKB vs MDR data collection
for i = 1:n
    NFKB_SA = NFKB_param(i);
    for j = 1:m
        MDR = damrate(j);
        ERiQ_SA
        life(i,j) = duration;
    end
end

%% Plot
X = 1.8E-3:0.1E-3:2.6E-3; % x-axis value - MDR
Y = {'-80%' '-60%' '-40%' '-20%' '0%' '20%' '40%' '60%' '80%'}; % y-axis value - NFKB

figure(1)
m = bar3(life);
ylabel('% Change in Rel. NFKB')
xlabel('Rel. Mito. Damage')
zlabel('Model Lifetime')
title('NFkB vs. Mitochondrial Damage')
set(gca,'XTickLabel',X) 
set(gca,'YTickLabel',Y) 
for k = 1:length(m)
    zdata = m(k).ZData;
    m(k).CData = zdata;
    m(k).FaceColor = 'interp';
end
hcb1 = colorbar;
pbaspect([1 1 0.3])
