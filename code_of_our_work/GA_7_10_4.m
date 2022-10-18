clear all
clc
tic

lb=[200,  0,  0,  0,  0,  0,200,  0,  0,  0,  0,  0];
ub=[900,450,450,950,550,550,900,450,450,950,550,550];
IntCon=[1,2,3,4,5,6,7,8,9,10,11,12];
nonlcon=[];
Ldata=[0:5:170,171:200];
para_result=zeros(12,length(Ldata));
para_resultbest=zeros(12,length(Ldata));
Rfin751=zeros(1,length(Ldata));
d=zeros(12,1);%实验参数赋值
d(1)=10^12;
d(2)=1.1;
d(3)=10^(-7);
d(4)=0.025;%分开xz
d(5)=0.025;
d(8)=10^(-15);
d(9)=10^(-10);
d(10)=10^(-10);
d(11)=10^(-10);
d(12)=10^(-10);
A=[0,0,0,1,1,1,0,0,0,0,0,0;
    -1,1,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,1,1,1;
    0,0,0,0,0,0,-1,1,0,0,0,0];
b=[999;0;999;0];
for j=1:length(Ldata)
    d(6)=0.4*10^(-0.2*(Ldata(j))/20);
    d(7)=0.4*10^(-0.2*(Ldata(j))/20);
    tot=25;
    fun=@(x)mdi(x,d);
    for i=1:tot
        [para_result(:,j)]=ga(fun,12,A,b,[],[],lb,ub,nonlcon,IntCon);
        mmp=-fun(para_result(:,j))/d(1);
        if mmp>Rfin751(j)
            Rfin751(j)=mmp;
            para_resultbest(:,j)=para_result(:,j);
        end
    end
end
save('/share/home/hlyin/bjl/data_final/para_resultbest_7.mat','para_resultbest')
save('/share/home/hlyin/bjl/data_final/Ldata_7.mat','Ldata')
save('/share/home/hlyin/bjl/data_final/Rfin75_7.mat','Rfin751')
toc