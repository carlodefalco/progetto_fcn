clear all
clc
nin=15;
nd=0;
nplot=50;

%out = A96_test_data (nin, nd, nplot)

x=linspace(-2,2,50)';
% y=x+2*x.^2.*sin(1./x);
% y=sin(x);
y=x.*sin(1./x);
Q=[x,y];
xd=zeros(nin,nd+1);
yd=zeros(nin,nd+1);
for i=2:nd+1
    xd(:,i-1)=out{1,i}(1,:);
    yd(:,i-1)=out{1,i}(2,:);
end
D=[xd,yd];

r=numel(Q(:,1))-1;
Wq=ones(1,r+1); %vettore valido per 
Wq(1,25)=0;
Wq(1,r+1)=1;
I=[0 1 3 4 5 6];
n=10;
p=3;

if nd==0
    s=-1
else
    s=nd;
end
Wd=ones(s+1,1);

[U,P]=WCLeastSquaresCurve(Q,r,Wq,D,s,I,Wd,n,p);

   figure
hold on, box on
plot(Q(:,1),Q(:,2),'o')

   [CC]=Grafico(n,p,U,P);
plot(P(:,1),P(:,2),'--')
title('x*sin(1/x) with f(x=0) constrained','fontsize',20)
legend('f(x)','Algorithm','location','southeast')
xlabel('x','Fontsize',24)
axis([-2 2 -0.4 1.2])
ylabel('f(x)','Fontsize',24)