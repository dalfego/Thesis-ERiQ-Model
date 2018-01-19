function MD = f_MD(MFUNCT,ROS)

global MDR MDR_SA

ScaleFactor = .8;

%MD = MDR_SA*((abs(MFUNCT+(ROS)))*MDR);
MD = MDR_SA*(((abs(MFUNCT+ROS))*MDR)+((ROS-ScaleFactor)*0.0001));


end