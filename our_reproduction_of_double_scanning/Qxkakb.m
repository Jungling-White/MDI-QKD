function Qxkakb=Qxkakb(pd,edx,etaa,etab,ka,kb)

Qxkakb=[0,0];

if (ka~=0)&&(kb~=0)
    I0 = besseli(0,sqrt(ka*etaa*kb*etab));
    I01 = besseli(0,sqrt(ka*etaa*kb*etab)/2);
else
    I0 =1;
    I01=1;
end

y=(1-pd)*exp(-(ka*etaa+kb*etab)/4);
Qc=2*y^2*(y^2-2*y*I01+I0);
Qe=2*y^2*(1+y^2-2*y*I01);
Qxkakb(1)=Qc+Qe;
Qxkakb(2)=(edx*Qc+(1-edx)*Qe);
%Qxkakb(2)等于Exkakb乘Qxkakb

end
