function [Cz,Az] = f_ROS_feedback(ROS,Cz,Uz)

gz = .01;

Cz = -ROS-Cz;
Az = gz*Uz;

end