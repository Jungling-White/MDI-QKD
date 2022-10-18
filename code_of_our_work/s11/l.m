function l=l(para,N,f,pd,edz,edx,etaa,etab,ecor,esec,e1,ec,ePA)
%para有16个元素，光强u,v,w,0和对应的概率。

QEz=Qzkakb(para,pd,edz,etaa,etab,para(1),para(9));
Qzuu=QEz(1);
Ezuu=QEz(2)/QEz(1);
nzuu=N*para(5)*para(13)*Qzuu;
%mzuu=N*para(5)^2*Ezuu*Qzuu;
%Sz0u=fSz0u(para,N,pd,etaa,etab,edz,esec);
SLEz11=fSz11(para,N,pd,etaa,etab,edz,esec);

phiz11L=fphiz11(para,N,SLEz11,pd,etaa,etab,edx,esec);
phiz11=phiz11L(1);
Sz11=phiz11L(3);
C=log2(2/ecor)+2*log2(2/(e1*ec))+2*log2(1/(2*ePA));
H1=-phiz11*log2(phiz11)-(1-phiz11)*log2(1-phiz11);
H2=-Ezuu*log2(Ezuu)-(1-Ezuu)*log2(1-Ezuu);
l(1)=Sz11*(1-H1)-nzuu*f*H2-C;
l(2)=Sz11;
l(3)=phiz11;
if l(1)<0
    l(1)=0;
end
if imag(l(1))~=0
    l(1)=0;
end
if phiz11>0.5
    l(1)=0;
end
if Ezuu>0.5
    l(1)=0;
end
end