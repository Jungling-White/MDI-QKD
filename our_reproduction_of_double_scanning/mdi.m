function R=mdi(para,d)

para1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
for i=1:3
    para1(i)=para(i)/1000;
end
para1(4)=0;
for j=1:3
    para1(j+4)=para(j+3)/1000;
end
para1(8)=1-para1(5)-para1(6)-para1(7);

para1(9)=para(7)/1000;
para1(10)=(para(2)/(para(8)/1000))/1000;
para1(11)=(para(3)/(para(9)/1000))/1000;
% para1(10)=para(8)/1000;
% para1(11)=para(9)/1000;

para1(12)=0;
for n=1:3
    para1(n+12)=para(n+9)/1000;
end
para1(16)=1-para1(13)-para1(14)-para1(15);

R=-l(para1,d(1),d(2),d(3),d(4),d(5),d(6),...
    d(7),d(8),d(9),d(10),d(11),d(12));

end
