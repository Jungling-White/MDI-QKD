function Sz0mu=fSz0mu(para,pd,edz,etaa,etab,N,epsilon)

Qz0u=Qzkakb(pd,edz,etaa,etab,0,para(9));
nz0mu=N*para(8)*para(13)*Qz0u(1);
ELnz0mu=chernoff(nz0mu,1,0,epsilon);

ELnz=exp(-para(1))*para(5)/para(8)*ELnz0mu;
nz=chernoff(ELnz,1,1,epsilon);

Sz0mu=nz;

end
