clear all
clc
tic
%完全按照论文中的方法，没有随机选样过程
lb=[100,  1,  1,100,  1, 50,100, 500, 500,100,  1, 50];
ub=[800,400,150,900,300,400,800,3000,3000,900,300,400];
IntCon=[1,2,3,4,5,6,7,8,9,10,11,12];
nonlcon=[];
Ldata=[0:10:90,91:1:110];
para_result=zeros(12,length(Ldata));
para_resultbest=zeros(12,length(Ldata));
Rfin=zeros(1,length(Ldata));
d=zeros(10,1);%实验参数赋值
d(1)=10^10;
d(2)=1.1;
d(3)=10^(-7);
d(4)=0.015;%分开xz
d(5)=0.015;
d(6)=10^(-15);
d(7)=10^(-10)*26;
d(8)=10^(-10);
d(9)=10^(-10);
d(10)=10^(-10);
A=[0,0,0,1,1,1,0,0,0,0,0,0;
    0,-1,1,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,1,1,1;
    0,0,0,0,0,0,0,1,-1,0,0,0];%va/vb<wa/wb
b=[999;0;999;0];

p=parpool(48);
parfor j=1:length(Ldata)
    l=zeros(2,1);
    l(1)=0.4*10^(-0.2*(Ldata(j))/20);
    l(2)=0.4*10^(-0.2*(Ldata(j))/20);
    tot=3;
    fun=@(x)mdi1(x,d,l);
    for i=1:tot
        [para_result(:,j)]=ga(fun,12,A,b,[],[],lb,ub,nonlcon,IntCon);
        Mp=-fun(para_result(:,j))/d(1);
        if Mp>Rfin(j)
            Rfin(j)=Mp;
            para_resultbest(:,j)=para_result(:,j);
        end
    end
end
delete(p);

semilogy(Ldata,Rfin,'r')
save('/share/home/hlyin/bjl/data/Rfin_9_9_1_as.mat','Rfin')
save('/share/home/hlyin/bjl/data/Ldata_9_9_1_as.mat','Ldata')
save('/share/home/hlyin/bjl/data/para_resultbest_9_9_1_as.mat','para_resultbest')
toc

function R=mdi1(para,d,l)

d1(1)=d(1);
d1(2)=d(2);
d1(3)=d(3);
d1(4)=d(4);
d1(5)=d(5);
d1(6)=l(1);
d1(7)=l(2);
d1(8)=d(6);
d1(9)=d(7);
d1(10)=d(8);
d1(11)=d(9);
d1(12)=d(10);
R1=mdi(para,d1);
R=R1(1);

end
