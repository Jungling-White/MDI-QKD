function phiz11=fphiz11(para,ESx11,N,M,H,epsilon)

phiz11=[0,0,0];

Etx11=M-exp(-para(3)-para(11))*N*para(7)*para(15)*H/2;
ESz11=para(1)*para(9)*para(5)*para(13)/(para(3)*para(11)*para(7)*para(15))...
    *exp((para(3)-para(1))+(para(11)-para(9)))*ESx11;
Sz11=chernoff(ESz11,1,1,epsilon);
% tx11=chernoff(Etx11*ESz11/ESx11,0,1,epsilon);
% Sx11=chernoff(ESx11,1,1,epsilon);
% 
% phiz11(1)=tx11/Sz11;
%tx11=chernoff(Etx11,0,1,epsilon);
Sx11=chernoff(ESx11,1,1,epsilon);

phiz11(1)=chernoff(Sz11*(Etx11/ESx11),0,1,epsilon)/Sz11;
phiz11(2)=Sz11;
phiz11(3)=Sx11;
