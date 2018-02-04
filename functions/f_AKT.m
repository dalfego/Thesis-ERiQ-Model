function AKT = f_AKT(PTEN,ROS)

global AKT_SA

GF = 0.1;
AKT = AKT_SA*(GF+PTEN+(ROS/5));

end
