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

pten = x; akt = x; nfkb = x; p53 = x; ampk = x; pgc = x; mtor = x; 
auto = x; foxo = x; sirt = x; hif = x; menzy = x; glyenz = x;
atpg = x; mdam = x; life = x;
pten2 = x; akt2 = x; nfkb2 = x; p532 = x; ampk2 = x; pgc2 = x; mtor2 = x; 
auto2 = x; foxo2 = x; sirt2 = x; hif2 = x; menzy2 = x; glyenz2 = x;
atpg2 = x; mdam2 = x; life2 = x;

for i = 1:n
    ERiQ_GSA_variation_drug
    pten(i) = PTEN_SA;        
    akt(i) = AKT_SA;      
    nfkb(i) = NFKB_SA;       
    p53(i) = P53_SA;          
    ampk(i) = AMPK_SA;       
    pgc(i) = PGC1a_SA;      
    mtor(i) = MTOR_SA;        
    auto(i) = AUTO_SA;       
    foxo(i) = FOXO_SA;        
    sirt(i) = SIRT_SA;       
    hif(i) = HIF_SA;          
    menzy(i) = MENZY_SA;      
    glyenz(i) = GLYENZ_SA; 
    ERiQ_GSA_variation_drug2
    atpg(i) = ATPg2(end);      
    mdam(i) = MDAMAGE2(end);   
    life(i) = duration2;       
    colNames = {'PTEN','AKT','NFKB','P53','AMPK','PGC1a','MTOR', ...
    'Autophagy','FOXO','SIRT','HIF','MENZY','GLYENZ','ATPg', ...
    'MDAMAGE','LIFETIME'};
    %GSA is for actual SA parameter scaling factor
    GSA_array = [pten akt nfkb p53 ampk pgc mtor auto foxo sirt hif menzy glyenz atpg mdam life];
    GSA = array2table(GSA_array,'VariableNames',colNames);
end
