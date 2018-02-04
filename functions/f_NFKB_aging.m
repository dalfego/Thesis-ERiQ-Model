function NFKB = f_NFKB_aging(AKT,ROS,MTOR,t)

global NFKB_SA

NFKB = NFKB_SA*((AKT+(ROS*0.25)+(MTOR*0.25))-(t*.0005));

end
