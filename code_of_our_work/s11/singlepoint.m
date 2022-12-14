% Ldata=80;
% para_result=[634,41,129,4,372,419,711,201,127,296,437,170];


d=zeros(12,1);%实验参数赋值
d(1)=10^10;
d(2)=1.1;
d(3)=10^(-7);
d(4)=0.015;%分开xz
d(5)=0.015;
d(8)=10^(-15);
d(9)=10^(-10);
d(10)=10^(-10);
d(11)=10^(-10);
d(12)=10^(-10);
R=zeros(3,length(Ldata));

for j=1:length(Ldata)
    d(6)=0.4*10^(-0.2*(Ldata(j))/20);
    d(7)=0.4*10^(-0.2*(Ldata(j))/20);

    fun=@(x)mdi(x,d);

    R(:,j)=-fun(para_result(:,j));
end
