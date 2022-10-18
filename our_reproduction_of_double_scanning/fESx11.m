function ESx11=fESx11(para,pd,edx,etaa,etab,N,M,H,epsilon)

A=N*para(7)*para(15)*para(3)*para(11)*exp(-para(3)-para(11))/...
    (para(2)*para(10)*para(3)*para(11)*(para(2)-para(3)));
a=para(2)^2*para(10)*exp(para(3)+para(11))/(N*para(7)*para(15));
b=para(3)^2*para(11)*exp(para(10))/(N*para(14)*para(8));
c=para(3)^2*para(11)*exp(para(2))/(N*para(6)*para(16));
d=para(3)^2*para(11)*exp(para(2)+para(10))/(N*para(6)*para(14));
f=para(3)^2*para(11)/(N*para(8)*para(16));
m=para(2)^2*para(10)*exp(para(3)+para(11))/(N*para(7)*para(15));
h=para(2)^2*para(10);
%A为总体因子，a为EUn~xww的系数，b为EUnxv0的系数，c为EUnx0v的系数；
%d为ELnxvv的系数，f为ELnx00的系数，m为M的系数，h为H的系数。

MQxww=Qxkakb(pd,edx,etaa,etab,para(3),para(11));
nxww=N*para(7)*para(15)*MQxww(1);
mxww=N*para(7)*para(15)*MQxww(2);
nCxww=nxww-mxww;
%第一个变量

MQxv0=Qxkakb(pd,edx,etaa,etab,para(2),0);
nxv0=N*para(6)*para(16)*MQxv0(1);
%第3个变量

MQx0v=Qxkakb(pd,edx,etaa,etab,0,para(10));
nx0v=N*para(8)*para(14)*MQx0v(1);
%第2个变量

MQxvv=Qxkakb(pd,edx,etaa,etab,para(2),para(10));
nxvv=N*para(6)*para(14)*MQxvv(1);
%第四个变量

MQx00=Qxkakb(pd,edx,etaa,etab,0,0);
nx00=N*para(8)*para(16)*MQx00(1);
%第五个变量

if a<b
    b1=a;
    c1=nCxww;
    b2=b;
    c2=nx0v;
else
    b1=b;
    c1=nx0v;
    b2=a;
    c2=nCxww;
end
if b2<c
    b3=c;
    c3=nxv0;
else
    if b1<c
        b3=b2;
        c3=c2;
        b2=c;
        c2=nxv0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=c;
        c1=nxv0;
    end
end
FL=b1*chernoff(c1+c2+c3,1,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,1,0,epsilon)+...
    (b3-b2)*chernoff(c3,1,0,epsilon);
%联合约束
a6=0;
if d<f
    b1=d;
    c1=nxvv;
    b2=f;
    c2=nx00;
else
    b1=f;
    c1=nx00;
    b2=d;
    c2=nxvv;
end
if b2<a6
    b3=a6;
    c3=0;
else
    if b1<a6
        b3=b2;
        c3=c2;
        b2=a6;
        c2=0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a6;
        c1=0;
    end
end
FU=b1*chernoff(c1+c2+c3,0,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,0,0,epsilon)+...
    (b3-b2)*chernoff(c3,0,0,epsilon);
%前三个变量采用联合约束，后两个变量补一个零系数乘零变量继续使用联合约束

ESx11=A*(FL-FU+m*M-h*H);

end
