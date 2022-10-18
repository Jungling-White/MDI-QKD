function phiz11=fphiz11(para,N,SLEz11,pd,etaa,etab,edx,esec)

phiz11=[0,0,0];

epsilon=esec/21;%失败概率相关计数

% YLEx11=SLEz11/(N*(para(1)^2*para(5)^2*exp(-2*para(1))+...
%     para(2)^2*para(6)^2*exp(-2*para(2))+...
%     2*para(1)*para(2)*exp(-para(1)-para(2))*para(5)*para(6)));
% SLEx11=N*para(7)^2*para(3)^2*exp(-2*para(3))*YLEx11;
% Sx11=chernoff(SLEx11,1,1,epsilon);
% Sz11=chernoff(SLEx11+SLEz11,1,1,epsilon)-Sx11;
% 
% a2=exp(-2*para(3))*para(7)^2/para(8)^2;
% a3=exp(-para(3))*para(7)/para(8);
% 
% QExww=Qxkakb(para,pd,edx,etaa,etab,para(3),para(3));
% x1=N*para(7)^2*QExww(2);
% QEx00=Qxkakb(para,pd,edx,etaa,etab,para(4),para(4));
% y2=chernoff(N*para(8)^2*QEx00(2),0,0,epsilon);
% x2=chernoff(a2*y2,0,1,epsilon);
% QExw0=Qxkakb(para,pd,edx,etaa,etab,para(3),para(4));
% y3=N*para(7)*para(8)*QExw0(2);
% QEx0w=Qxkakb(para,pd,edx,etaa,etab,para(4),para(3));
% y4=N*para(8)*para(7)*QEx0w(2);
% x3=chernoff(a3*chernoff(y3+y4,1,0,epsilon),1,1,epsilon);
% tx11=x1+x2-x3;
% 
% phiz11(1)=tx11/Sz11+randomsample(Sz11,Sx11,tx11/Sz11,epsilon);

YLEx11=SLEz11/(N*(para(1)*para(9)*para(5)*para(13)*exp(-para(1)-para(9))+...
    para(2)*para(10)*para(6)*para(14)*exp(-para(2)-para(10))+...
    para(1)*para(2)*exp(-para(1)-para(2))*para(5)*para(6)+...
    para(9)*para(10)*exp(-para(9)-para(10))*para(13)*para(14)));
SLEx11=N*para(7)*para(15)*para(3)*para(11)*exp(-para(3)-para(11))*YLEx11;
Sx11=chernoff(SLEx11,1,1,epsilon);
Sz11=chernoff(SLEx11+SLEz11,1,1,epsilon)-Sx11;

a1=1;
a2=1/2*exp(-para(3)-para(11))*para(7)*para(15)/(para(8)*para(16));
a3=1/2*exp(-para(3))*para(7)/para(8);
a4=1/2*exp(-para(11))*para(15)/para(16);

QExww=Qxkakb(para,pd,edx,etaa,etab,para(3),para(11));
x1=chernoff(N*para(7)*para(15)*QExww(2),0,0,epsilon);
QEx00=Qxkakb(para,pd,edx,etaa,etab,para(4),para(12));
x2=chernoff(N*para(8)*para(16)*QEx00(1),0,0,epsilon);
QEx0w=Qxkakb(para,pd,edx,etaa,etab,para(4),para(11));
x3=chernoff(N*para(8)*para(15)*QEx0w(1),1,0,epsilon);
QExw0=Qxkakb(para,pd,edx,etaa,etab,para(3),para(12));
x4=chernoff(N*para(7)*para(16)*QExw0(1),1,0,epsilon);

tx11=chernoff(SLEz11*(a1*x1+a2*x2-a3*x3-a4*x4)/SLEx11,...
    0,1,epsilon);

phiz11(1)=tx11/Sz11;

phiz11(2)=Sx11;
phiz11(3)=Sz11;

end