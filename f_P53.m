function [P53s,P53a,P53] = f_P53(AKT,NFKB,ROS,Ax)

global P53_Base P53_Act P53_SA

%P53Act = 1 in script {< 1 = no activation, modify this between 0.1 and 4};

P53s = 0.3*(P53_Base-(AKT)-(NFKB)+(ROS*0.5))*P53_Act;
P53a = Ax-P53s;             %{y, output tracks r}
P53 = P53_SA*(P53s + P53a);

end