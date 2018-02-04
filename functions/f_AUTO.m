function AUTOPHAGY = f_AUTO(MTOR,FOXO,ROS,P53)

global AUTO_SA
Scale1 = 0.001;

AUTOPHAGY = AUTO_SA*(Scale1*((1./MTOR)+(FOXO*0.5)+ROS+P53));

end
