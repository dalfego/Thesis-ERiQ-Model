function NFKB = f_NFKB(AKT,ROS,MTOR)

global NFKB_SA

NFKB = NFKB_SA*(AKT+(ROS*0.25)+(MTOR*0.25))*1;

end
