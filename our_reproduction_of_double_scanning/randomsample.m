function randomsample=randomsample(n,k,lambda,epsilon)

g=2*(n+k)/(n*k)*log(sqrt(n+k)/sqrt(2*pi*n*k*lambda*(1-lambda))...
    *epsilon);

if n>=k
    A=(1-2*lambda)*n*g/(n+k);
    B=sqrt((n^2*g^2)/(n+k)^2+4*lambda*g-4*lambda^2*g);
    C=2+2*(n^2*g)/(n+k)^2;
    randomsample=(A+B)/C; 
else
    A=(1-2*lambda)*k*g/(n+k);
    B=sqrt((k^2*g^2)/(n+k)^2+4*lambda*g-4*lambda^2*g);
    C=2+2*(k^2*g)/(n+k)^2;
    randomsample=(A+B)/C; 
end

end