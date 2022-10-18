function Sz0u=fSz0u(para,N,pd,etaa,etab,ed,esec)

QEz=Qzkakb(para,pd,ed,etaa,etab,0,para(9));
Qz0u=QEz(1);
nz0u=N*para(8)*para(13)*Qz0u;
epsilon=esec/21;%失败概率相关计数

nLEz0u=chernoff(nz0u,1,0,epsilon);
nLEz=(1+exp(-para(1))*para(5)/para(8))*nLEz0u;
nLz=chernoff(nLEz,1,1,epsilon);%观测值
Sz0u=nLz-nz0u;

end