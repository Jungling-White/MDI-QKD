function l=l(para,N,f,pd,edz,edx,etaa,etab,ecor,esec,e1,ec,ePA)
%para有16个元素，光强u,v,w,0和对应的概率。

l=[0,0,0,0];
%l(1)码长;l(2)相位错误;l(3)Sz11观测值;l(4)Sx11观测值;
epsilon=esec/26;%失败概率相关计数

MQxww=Qxkakb(pd,edx,etaa,etab,para(3),para(11));
mxww=N*para(7)*para(15)*MQxww(2);
ML=chernoff(mxww,1,0,epsilon);
MU=chernoff(mxww,0,0,epsilon);
%双扫描方法涉及的第一个变量M期望值的上下界

MQxw0=Qxkakb(pd,edx,etaa,etab,para(3),0);
MQx0w=Qxkakb(pd,edx,etaa,etab,0,para(11));
MQx00=Qxkakb(pd,edx,etaa,etab,0,0);
nxw0=N*para(7)*para(16)*MQxw0(1);
nx0w=N*para(8)*para(15)*MQx0w(1);
nx00=N*para(8)*para(16)*MQx00(1);
a1=exp(para(3))/(N*para(7)*para(16));
a2=exp(para(11))/(N*para(8)*para(15));
a3=0;
if a1<a2
    b1=a1;
    c1=nxw0;
    b2=a2;
    c2=nx0w;
else
    b1=a2;
    c1=nx0w;
    b2=a1;
    c2=nxw0;
end
if b2<a3
    b3=a3;
    c3=0;
else
    if b1<a3
        b3=b2;
        c3=c2;
        b2=a3;
        c2=0;
    else
        b3=b2;
        c3=c2;
        b2=b1;
        c2=c1;
        b1=a3;
        c1=0;
    end
end
FL=b1*chernoff(c1+c2+c3,1,0,epsilon)+...
    (b2-b1)*chernoff(c2+c3,1,0,epsilon)+...
    (b3-b2)*chernoff(c3,1,0,epsilon);
HL=FL-chernoff(nx00,0,0,epsilon)/(N*para(8)*para(16));

a4=exp(para(3))/(N*para(7)*para(16));
a5=exp(para(11))/(N*para(8)*para(15));
a6=0;
if a4<a5
    b1=a4;
    c1=nxw0;
    b2=a5;
    c2=nx0w;
else
    b1=a5;
    c1=nx0w;
    b2=a4;
    c2=nxw0;
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
HU=FU-chernoff(nx00,1,0,epsilon)/(N*para(8)*para(16));
%双扫描方法涉及的第二个变量H的上下界，上下界由联合约束得到
if MU<ML
    l(1)=0;
    l(2)=0;
    l(3)=0;
    l(4)=0;
    return
end
if HU<HL
    l(1)=0;
    l(2)=0;
    l(3)=0;
    l(4)=0;
    return
end

% %original method begin
% para_result1=zeros(2,1);
% rest(1:16)=para;
% rest(9+8)=N;
% rest(10+8)=f;
% rest(11+8)=pd;
% rest(12+8)=edz;
% rest(13+8)=edx;
% rest(14+8)=etaa;
% rest(15+8)=etab;
% rest(16+8)=ecor;
% rest(17+8)=esec;
% rest(18+8)=e1;
% rest(19+8)=ec;
% rest(20+8)=ePA;
% rest(21+8)=ML;
% rest(22+8)=MU;
% rest(23+8)=HL;
% rest(24+8)=HU;
% 
% 
% myfun=@(x)l1(x,rest);
% Rfin=myfun(para_result1(:,1));
% 
% ESx11=fESx11(para,pd,edx,etaa,etab,N,ML,...
%     HU,epsilon);
% 
% phiz11z=fphiz11(para,ESx11,N,MU,...
%     HL,epsilon);
% %original method end

%begin double-scanning
ML1=0;
HL1=0;
MU1=300;
HU1=300;
lb=[ML1,HL1];
ub=[MU1,HU1];
IntCon=[1,2];
nonlcon=[];
para_result=zeros(2,1);
para_result1=zeros(2,1);
para_resultbest1=zeros(2,1);
rest(1:16)=para;
rest(9+8)=N;
rest(10+8)=f;
rest(11+8)=pd;
rest(12+8)=edz;
rest(13+8)=edx;
rest(14+8)=etaa;
rest(15+8)=etab;
rest(16+8)=ecor;
rest(17+8)=esec;
rest(18+8)=e1;
rest(19+8)=ec;
rest(20+8)=ePA;
rest(21+8)=ML;
rest(22+8)=MU;
rest(23+8)=HL;
rest(24+8)=HU;

min=N;
tot=4;
for i=1:tot
    myfun=@(x)l1(x,rest);
    [para_result(:,1)]=ga(myfun,2,[],[],[],[],lb,ub,...
        nonlcon,IntCon);
    if min>=myfun(para_result(:,1))
        min=myfun(para_result(:,1));
        para_result1(:,1)=para_result(:,1);
    end
end
Rfin=myfun(para_result1(:,1));

para_resultbest1(1,1)=para_result1(1,1)*(MU-ML)/300+ML;
para_resultbest1(2,1)=para_result1(2,1)*(HU-HL)/300+HL;

ESx11=fESx11(para,pd,edx,etaa,etab,N,para_resultbest1(1,1),...
    para_resultbest1(2,1),epsilon);

phiz11z=fphiz11(para,ESx11,N,para_resultbest1(1,1),...
    para_resultbest1(2,1),epsilon);
%end double-scanning
phiz11=phiz11z(1);
Sz11=phiz11z(2);
Sx11=phiz11z(3);

l(1)=Rfin;
l(2)=phiz11;
l(3)=Sz11;
l(4)=Sx11;
%searching condition
% if l(1)==N
%     l(1)=0;
% end
if l(1)<0
    l(1)=0;
end
if imag(l(1))~=0
    l(1)=0;
end
if phiz11>0.5
    l(1)=0;
end

end
