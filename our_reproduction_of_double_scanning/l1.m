function l1=l1(x,rest)
%x为M与H，
%rest为para(为光强与概率参数)，N,f,pd,edz,edx,
%etaa,etab,ecor,esec,e1,ec,ePA。

para=rest(1:16);
N=rest(9+8);
f=rest(10+8);
pd=rest(11+8);
edz=rest(12+8);
edx=rest(13+8);
etaa=rest(14+8);
etab=rest(15+8);
ecor=rest(16+8);
esec=rest(17+8);
e1=rest(18+8);
ec=rest(19+8);
ePA=rest(20+8);
% %original method begin
% ML=rest(21+8);
% MU=rest(22+8);
% HL=rest(23+8);
% HU=rest(24+8);
% 
% epsilon=esec/26;%失败概率相关计数
% QEz=Qzkakb(pd,edz,etaa,etab,para(1),para(9));
% Qzuu=QEz(1);
% Ezuu=QEz(2)/QEz(1);
% nzuu=N*para(5)*para(13)*Qzuu;
% 
% Sz0u=fSz0mu(para,pd,edz,etaa,etab,N,epsilon);
% ESx11=fESx11(para,pd,edx,etaa,etab,N,ML,HU,epsilon);
% 
% phiz11z=fphiz11(para,ESx11,N,MU,HL,epsilon);
% phiz11=phiz11z(1);
% Sz11=phiz11z(2);
% C=log2(8/ecor)+2*log2(2/(e1*ec))+2*log2(1/(2*ePA));
% H1=-phiz11*log2(phiz11)-(1-phiz11)*log2(1-phiz11);
% H2=-Ezuu*log2(Ezuu)-(1-Ezuu)*log2(1-Ezuu);
% %original method end

%double scanning method begin
M=x(1)*(rest(22+8)-rest(21+8))/300+rest(21+8);
H=x(2)*(rest(24+8)-rest(23+8))/300+rest(23+8);

epsilon=esec/26;%失败概率相关计数
QEz=Qzkakb(pd,edz,etaa,etab,para(1),para(9));
Qzuu=QEz(1);
Ezuu=QEz(2)/QEz(1);
nzuu=N*para(5)*para(13)*Qzuu;

Sz0u=fSz0mu(para,pd,edz,etaa,etab,N,epsilon);
ESx11=fESx11(para,pd,edx,etaa,etab,N,M,H,epsilon);

phiz11z=fphiz11(para,ESx11,N,M,H,epsilon);
phiz11=phiz11z(1);
Sz11=phiz11z(2);
C=log2(8/ecor)+2*log2(2/(e1*ec))+2*log2(1/(2*ePA));
H1=-phiz11*log2(phiz11)-(1-phiz11)*log2(1-phiz11);
H2=-Ezuu*log2(Ezuu)-(1-Ezuu)*log2(1-Ezuu);
%double scanning method end

l1=(Sz0u+Sz11*(1-H1)-nzuu*f*H2-C);
if Ezuu>1/2
    l1=0;
end
%searching condition
% if l1<0
%     l1=N;
% end
% if imag(l1)~=0
%     l1=N;
% end
% if phiz11>0.5
%     l1=N;
% end

end
