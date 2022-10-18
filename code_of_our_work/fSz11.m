function Sz11=fSz11(para,N,pd,etaa,etab,edz,esec)

epsilon=esec/21;%失败概率相关计数
% A=para(1)*(para(1)^2*para(5)^2*exp(-2*para(1))+...
%     para(2)^2*para(6)^2*exp(-2*para(2))+...
%     2*para(1)*para(2)*exp(-para(1)-para(2))*para(5)*para(6))...
%     /(2*(para(1)*para(2)-para(2)^2));
% a=para(1)*exp(para(2))/(para(1)*para(2)-para(2)^2)*2*exp(para(2))/para(6)^2;
% b=para(2)/para(1)^2*exp(para(1))/(para(5)*para(8));
% c=2*(para(1)^2-para(2)^2)/(para(1)^2*para(2))/para(8)^2;
% d=(para(2)/para(1)^2*exp(para(1)+para(2))+para(1)*exp(para(2))/(para(1)*para(2)...
%     -para(2)^2)*para(2)^2/para(1)^2*exp(para(1)))/(para(5)*para(6));
% e=((para(1)^2-para(2)^2)/(para(1)^2*para(2))*exp(para(2))+para(1)*exp(para(2))...
%     /(para(1)*para(2)-para(2)^2)*(para(1)^2-para(2)^2)/para(1)^2)/(para(6)*para(8));
% a1=A*a;
% a2=A*b;
% a3=A*c;
% a4=A*d;
% a5=A*e;

Out=(para(1)*para(9)*para(5)*para(13)*exp(-para(1)-para(9))+...
    para(2)*para(10)*para(6)*para(14)*exp(-para(2)-para(10))+...
    para(1)*para(2)*exp(-para(1)-para(2))*para(5)*para(6)+...
    para(9)*para(10)*exp(-para(9)-para(10))*para(13)*para(14))/2;
A=para(9)/(para(9)*para(10)-para(10)^2);
B=para(1)/(para(1)*para(2)-para(2)^2);
C=para(2)^2/para(1)^2;
D=para(10)^2/para(9)^2;
E=(para(1)^2-para(2)^2)/para(1)^2;
F=(para(9)^2-para(10)^2)/para(9)^2;
A1=A*B*exp(para(2)+para(10))/(para(6)*para(14));
A2=A1;
B1=A*B*C*exp(para(1)+para(10))/(para(5)*para(14));
B2=B*A*D*exp(para(2)+para(9))/(para(6)*para(13));
C1=A*B*E*exp(para(10))/(para(8)*para(14));
C2=B*A*F*exp(para(2))/(para(6)*para(16));
D1=A*D*exp(para(2)+para(9))/para(2)/(para(6)*para(13));
D2=B*C*exp(para(1)+para(10))/para(10)/(para(5)*para(14));
E1=A*D*exp(para(9))/para(2)/(para(8)*para(13));
E2=B*C*exp(para(1))/para(10)/(para(5)*para(16));
F1=A*F*exp(para(2))/para(2)/(para(6)*para(16));
F2=B*E*exp(para(10))/para(10)/(para(8)*para(14));
G1=A*F/para(2)/(para(8)*para(16));
G2=B*E/para(10)/(para(8)*para(16));
a1=A1+A2;
a2=G1+G2;
a3=E1;
a4=E2;
a5=B1+D2;
a6=B2+D1;
a7=C1+F2;
a8=C2+F1;
%系数

QEzvv=Qzkakb(para,pd,edz,etaa,etab,para(2),para(10));
Qzvv=QEzvv(1);
nzvv=N*para(6)*para(14)*Qzvv;
%第一个变量

QEz00=Qzkakb(para,pd,edz,etaa,etab,0,0);
Qz00=QEz00(1);
nz00=N*para(8)*para(16)*Qz00;
%第二个变量

QEz0u=Qzkakb(para,pd,edz,etaa,etab,0,para(9));
Qz0u=QEz0u(1);
nz0u=N*para(8)*para(13)*Qz0u;
%第三个变量

QEzu0=Qzkakb(para,pd,edz,etaa,etab,para(1),0);
Qzu0=QEzu0(1);
nzu0=N*para(5)*para(16)*Qzu0;
%第四个变量

QEzuv=Qzkakb(para,pd,edz,etaa,etab,para(1),para(10));
Qzuv=QEzuv(1);
nzuv=N*para(5)*para(14)*Qzuv;
%第5个变量

QEzvu=Qzkakb(para,pd,edz,etaa,etab,para(2),para(9));
Qzvu=QEzvu(1);
nzvu=N*para(6)*para(13)*Qzvu;
%第6个变量

QEz0v=Qzkakb(para,pd,edz,etaa,etab,0,para(10));
Qz0v=QEz0v(1);
nz0v=N*para(8)*para(14)*Qz0v;
%第7个变量

QEzv0=Qzkakb(para,pd,edz,etaa,etab,para(2),0);
Qzv0=QEzv0(1);
nzv0=N*para(6)*para(16)*Qzv0;
%第8个变量

if a1<a2
    b1=a1;
    c1=nzvv;
    b2=a2;
    c2=nz00;
else
    b1=a2;
    c1=nz00;
    b2=a1;
    c2=nzvv;
end
if b2<a3
    b3=a3;
    c3=nz0u;
else
    if b1<a3
        b3=b2;
        c3=c2;
        b2=a3;
        c2=nz0u;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a3;
        c1=nz0u;
    end
end
if b3<a4
    b4=a4;
    c4=nzu0;
else
    if b2<a4
        b4=b3;
        c4=c3;
        b3=a4;
        c3=nzu0;
    elseif b1<a4
        b4=b3;
        c4=c3;
        b3=b2;
        c3=c2;
        b2=a4;
        c2=nzu0;
    else
        b4=b3;
        c4=c3;
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a4;
        c1=nzu0;
    end
end
FL=b1*chernoff(c1+c2+c3+c4,1,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3+c4,1,0,epsilon)+...
    (b3-b2)*chernoff(c3+c4,1,0,epsilon)+(b4-b3)*chernoff(c4,1,0,epsilon);

if a5<a6
    b1=a5;
    c1=nzuv;
    b2=a6;
    c2=nzvu;
else
    b1=a6;
    c1=nzvu;
    b2=a5;
    c2=nzuv;
end
if b2<a7
    b3=a7;
    c3=nz0v;
else
    if b1<a7
        b3=b2;
        c3=c2;
        b2=a7;
        c2=nz0v;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a7;
        c1=nz0v;
    end
end
if b3<a8
    b4=a8;
    c4=nzv0;
else
    if b2<a8
        b4=b3;
        c4=c3;
        b3=a8;
        c3=nzv0;
    elseif b1<a8
        b4=b3;
        c4=c3;
        b3=b2;
        c3=c2;
        b2=a8;
        c2=nzv0;
    else
        b4=b3;
        c4=c3;
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a8;
        c1=nzv0;
    end
end
FU=b1*chernoff(c1+c2+c3+c4,0,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3+c4,0,0,epsilon)+...
    (b3-b2)*chernoff(c3+c4,0,0,epsilon)+(b4-b3)*chernoff(c4,0,0,epsilon);
%采用联合约束
Sz11=Out*(FL-FU);

end
