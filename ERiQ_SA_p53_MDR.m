%% Setup
global P53_Act NFKB_SA MDR ROS_SA2
ROS_SA2 = 1;
NFKB_SA = 1;

% Set range of desired parameters
p53_param = [0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0];
damrate = (1.5E-3):(0.1E-3):(2.4E-3);

n = numel(p53_param);
m = numel(damrate);

x = zeros(n,1);
p53 = x;
mdr = x;
life = x;

%% NFKB vs MDR data collection
for i = 1:n
    P53_Act = p53_param(i);
    for j = 1:m
        MDR = damrate(j);
        ERiQ_SA
        life(i,j) = duration;
    end
end

%% Plot
X1 = 1.5E-3:0.1E-3:2.4E-3; % x-axis value - MDR
Y1 = {'-60%' '-20%' '20%' '60%' '100%' '140%' '180%' '220%' '260%' '300%'}; % y-axis value - p53

figure(1)
m = bar3(life);
ylabel('Rel. p53')
xlabel('Rel. Mito. Damage')
zlabel('Model Lifetime')
title('p53 vs. Mitochondrial Damage')
set(gca,'XTickLabel', X1)
set(gca,'YTickLabel', Y1)
for k = 1:length(m)
    zdata = m(k).ZData;
    m(k).CData = zdata;
    m(k).FaceColor = 'interp';
end
hcb1 = colorbar;
pbaspect([1 1 0.3])
