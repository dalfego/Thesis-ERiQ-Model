function [Cx,Ax] = f_P53_feedback(P53a,Cx)

rx = 0;                       % { Zero P53a output if no change }
gx = 0.1;                     % { smooth close to 0 behaivior long term }
Ux = (rx+Cx);

Cx = -P53a-Cx;
Ax = gx*Ux;

end