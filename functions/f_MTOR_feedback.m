function [Cy,Ay] = f_MTOR_feedback(MTORa,Cy)

ry = 0;
gy = 0.1;
Uy = (ry+Cy);

Cy = -MTORa-Cy;
Ay = gy*Uy;

end
