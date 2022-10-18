function Qzkakb=Qzkakb(para,pd,edz,etaa,etab,ka,kb)

for i=1:4
   if para(i)==ka
       break
   end
end
for j=1:4
    if para(j)==kb
        break
    end
end

Qzkakb=[0,0];

I0 = besseli(0,sqrt(ka*etaa*kb*etab));
A=I0-(1-pd)*exp(-(ka*etaa+kb*etab)/2);
B=(1-(1-pd)*exp(-(ka*etaa)/2))*(1-(1-pd)*exp(-(kb*etab)/2));
Qc=2*4*(1-pd)^2*exp(-(ka*etaa+kb*etab)/2)*B;
Qe=2*4*pd*(1-pd)^2*exp(-(ka*etaa+kb*etab)/2)*A;
Qzkakb(1)=(Qc+Qe)/4;
Qzkakb(2)=(edz*Qc+(1-edz)*Qe)/4;
%Qzkakb(2)等于Ezkakb乘Qzkakb

end