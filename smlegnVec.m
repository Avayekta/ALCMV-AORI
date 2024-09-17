function [segnvector] = smlegnVec(eta_matrix)
%3x3 symmetrix matrix
a=eta_matrix;
alpha=a(1,1)+a(2,2)+a(3,3);
beta=(a(1,2)^2)+(a(1,3)^2)+(a(2,3)^2)-(a(1,1)*a(2,2))-(a(2,2)*a(3,3))-(a(3,3)*a(1,1));
gamma=a(1,1)*a(2,2)*a(3,3)+2*a(1,2)*a(2,3)*a(1,3)-a(1,1)*(a(2,3)^2)-(a(1,2)^2)*a(3,3)-(a(1,3)^2)*a(2,2);
p=-((3*beta)+(alpha^2))/3;
q=-(gamma+(2*(alpha^3)/27)+(alpha*beta)/3);
temp1=2*sqrt((abs(p)^3)/27);
phi=acos(-1*q*(1/temp1));
landa1=alpha/3+2*cos(phi/3)*sqrt(abs(p)/3);
landa2=alpha/3-2*cos((phi-pi)/3)*sqrt(abs(p)/3);
landa3=alpha/3-2*cos((phi+pi)/3)*sqrt(abs(p)/3);
landamin=min(abs(landa1),abs(landa2));
landamin=min(landamin,abs(landa3));
b1(1,1)=a(1,1)-landamin;
b1(2,2)=a(2,2)-landamin;
b1(3,3)=a(3,3)-landamin;
sbc1=(b1(1,1)*a(2,3)-a(1,3)*a(1,2))*a(1,3);
sbcs1=(a(1,2)^2-b1(1,1)*b1(2,2))*a(1,3);
if sbc1 ~= 0 | sbc1 ~= 0 
Q1=(b1(1,1)*a(2,3)-a(1,3)*a(1,2))/(a(1,2)^2-b1(1,1)*b1(2,2));
R1=1/Q1;
Pn1=-1*(a(2,3)*Q1+b1(3,3))/a(1,3);
den1=sqrt(Q1^2+Pn1^2+1);
segnvector=[Pn1;Q1;1]/den1;
else
sbc2=(b1(1,1)*b1(3,3)-a(1,3)^2)*a(1,2);
sbcs2=(a(1,2)*a(1,3)-b1(1,1)*a(2,3))*a(1,2);
if sbc2 ~= 0 | sbc2 ~= 0 
Q2=(b1(1,1)*b1(3,3)-a(1,3)^2)/(a(1,2)*a(1,3)-b1(1,1)*a(2,3));
R2=1/Q2;
Pn2=-1*(b1(2,2)*Q2+a(2,3))/a(1,2);
den2=sqrt(Q2^2+Pn2^2+1);
segnvector=[Pn2;Q2;1]/den2;
else
sbc3=(a(1,2)*b1(3,3)-a(2,3)*a(1,3))*b1(1,1);
sbcs3=(b1(2,2)*a(1,3)-a(1,2)*a(2,3))*b1(1,1);
if sbc3 ~= 0 | sbc3 ~= 0 
Q3=(a(1,2)*b1(3,3)-a(2,3)*a(1,3))/(b1(2,2)*a(1,3)-a(1,2)*a(2,3));
R3=1/Q3;
Pn3=-1*(a(1,2)*Q3+a(1,3))/b1(1,1);
den3=sqrt(Q3^2+Pn3^2+1);
segnvector=[Pn3;Q3;1]/den3;
else
sbc4=(a(1,2)*a(2,3)-a(1,3)*b1(2,2))*a(2,3);
sbcs4=(b1(1,1)*b1(2,2)-a(1,2)^2)*a(2,3);
if sbc4 ~= 0 | sbc4 ~= 0 
P4=(a(1,2)*a(2,3)-a(1,3)*b1(2,2))/(b1(1,1)*b1(2,2)-a(1,2)^2);
R4=1/P4;
Qn4=-1*(a(1,3)*P4+b1(3,3))/a(2,3);
den4=sqrt(Q4^2+P4^2+1);
segnvector=[P4;Q4;1]/den4;
else
sbc5=(a(1,2)*b1(3,3)-a(1,3)*a(2,3))*b1(2,2);
sbcs5=(b1(1,1)*a(2,3)-a(1,2)*a(1,3))*b1(2,2);
if sbc5 ~= 0 | sbc5 ~= 0 
P5=(a(1,2)*b1(3,3)-a(1,3)*a(2,3))/(b1(1,1)*a(2,3)-a(1,2)*a(1,3));
R5=1/P5;
Qn5=-1*(a(1,2)*P5+a(2,3))/b1(2,2);
den5=sqrt(Qn5^2+P5^2+1);
segnvector=[P5;Qn5;1]/den5;   
else
sbc6=(b1(2,2)*b1(3,3)-a(2,3)^2)*a(1,2);
sbcs6=(a(1,2)*a(2,3)-b1(2,2)*a(1,3))*a(1,2);
if sbc6 ~= 0 | sbc6 ~= 0 
P6=(b1(2,2)*b1(3,3)-a(2,3)^2)/(a(1,2)*a(2,3)-b1(2,2)*a(1,3));
R6=1/P6;
Qn6=-1*(b1(1,1)*P6+a(1,3))/a(1,2);
den6=sqrt(Qn6^2+P6^2+1);
segnvector=[P6;Qn6;1]/den6; 
else
sbc7=(a(1,3)*b1(2,2)-a(1,2)*a(2,3))*b1(3,3);
sbcs7=(b1(1,1)*a(2,3)-a(1,3))*b1(3,3);
if sbc7 ~= 0 | sbc7 ~= 0 
P7=(a(1,3)*b1(2,2)-a(1,2)*a(2,3))/(b1(1,1)*a(2,3)-a(1,3));
Q7=1/P7;
Rm7=-1*(a(1,3)*P7+a(2,3))/b1(3,3);
den7=sqrt(P7^2+Rm7^2+1);
segnvector=[P7;Q7;1]/den7;  
else
sbc8=(a(1,3)*a(2,3)-a(1,2)*b1(3,3))*a(2,3);
sbcs8=(b1(1,1)*b1(3,3)-a(1,3)^2)*a(2,3);
if sbc8 ~= 0 | sbc8 ~= 0 
P8=(a(1,3)*a(2,3)-a(1,2)*b1(3,3))/(b1(1,1)*b1(3,3)-a(1,3)^2);
Q8=1/P8;
Rm8=-1*(a(1,2)*P8+b1(2,2))/a(2,3);
den8=sqrt(Rm8^2+P8^2+1);
segnvector=[P8;1;Rm8]/den8;  
else
sbc9=(a(2,3)^2-b1(2,2)*b1(3,3))*a(1,3);
sbcs9=(a(1,2)*b1(3,3)-a(2,3)*a(1,3))*a(1,3);
if sbc9 ~= 0 | sbc9 ~= 0 
P9=(a(2,3)^2-b1(2,2)*b1(3,3))/(a(1,2)*b1(3,3)-a(2,3)*a(1,3));
Q9=1/P9;
Rm9=-1*(b1(1,1)*P9+a(1,2))/a(1,3);
den9=sqrt(Rm9^2+P9^2+1);
segnvector=[P9;1;Rm9]/den9;
end; end; end; end; end; end; end; end; end;
end


