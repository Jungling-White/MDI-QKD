function chernoff=chernoff(x,U,E,epsilon)
%U=0表示求上限，U=1表示求下限。E=0表示从观测值到期望值，E=1反之。

beta=log(1/epsilon);
triangleL=(beta+sqrt(8*beta*x+beta^2))/2;
triangleU=beta+sqrt(2*beta*x+beta^2);
deltaL=sqrt(2*beta/x);
deltaU=(beta+sqrt(8*beta*x+beta^2))/(2*x);
if E==0
    if U==0
        chernoff=x+triangleU;
    else
        chernoff=x-triangleL;
        if chernoff<0
            chernoff=0;
        end
    end
else
    if x==0
        chernoff=0;
        return
    end
    if U==0
        chernoff=x+deltaU*x;
    else
        chernoff=x-deltaL*x;
        if chernoff<0
            chernoff=0;
        end
    end
end
