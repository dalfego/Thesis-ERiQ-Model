% n = # of simulations
n = 10000;
x = zeros(n,1);
% original SA values = 1; this value is set to s in case changed
s = 1;

% REMOVE: ROS, MDR, GLU, PYR, NAD
% Make bar graphs to make sure there are no flipping points in lifetime
% increasing/decreasing parameters
% Increase: p53, Auto,
% Decrease: mTOR, AKT, NFkB
% What is the theoretical max life?
% GSA for standard model and increase for SIRT only

pten = x; akt = x; nfkb = x; p53 = x; mtor = x; auto = x;
atpg = x; mdam = x; life = x;
pten2 = x; akt2 = x; nfkb2 = x; p532 = x; mtor2 = x; auto2 = x;
atpg2 = x; mdam2 = x; life2 = x;

for i = 1:n
    ERiQ_Global_Analysis_OPTIMIZE
    pten(i) = PTEN_SA;        pten2(i) = 100*((PTEN_SA-s)/s);
    akt(i) = AKT_SA;          akt2(i) = 100*((AKT_SA-s)/s);
    nfkb(i) = NFKB_SA;        nfkb2(i) = 100*((NFKB_SA-s)/s);
    p53(i) = P53_SA;          p532(i) = 100*((P53_SA-s)/s); 
    mtor(i) = MTOR_SA;        mtor2(i) = 100*((MTOR_SA-s)/s);
    auto(i) = AUTO_SA;        auto2(i) = 100*((AUTO_SA-s)/s);
    atpg(i) = ATPg(end);      atpg2(i) = 100*((ATPg(end)-3.5358)/3.5358);
    mdam(i) = MDAMAGE(end);   mdam2(i) = 100*((MDAMAGE(end)-1.6559)/1.6559);
    life(i) = duration;       life2(i) = 100*((duration-752.9255)/752.9255);
    colNames = {'PTEN','AKT','NFKB','P53','MTOR','Autophagy','ATPg', ...
    'MDAMAGE','LIFETIME'};
    %GSA is for actual SA parameter scaling factor
    GSA_array = [pten akt nfkb p53 mtor auto atpg mdam life];
    GSA = array2table(GSA_array,'VariableNames',colNames);
    %GSA_percent is for the percent difference used in SA
    GSA_percent_array = [pten2 akt2 nfkb2 p532 mtor2 auto2 atpg2 mdam2 life2];
    GSA_percent = array2table(GSA_percent_array,'VariableNames',colNames);
end