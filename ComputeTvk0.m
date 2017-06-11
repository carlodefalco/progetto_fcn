function [Tv]=ComputeTvk0(n,m,Q,k)
for k=0:n
    alphak=zeros(n+1,m+1);
    q=zeros(n+4,m+1);
    for i=1:n
        q(i+2,1)=sqrt((Q(i+1,1,1)-Q(i,1,1))^2+(Q(i+1,1,2)-Q(i,1,2))^2+(Q(i+1,1,3)-Q(i,1,3))^2);
    end
    q(2,1)=2*q(3,1)-q(4,1);
    q(1,1)=2*q(2,1)-q(3,1);
    q(n+3,1)=2*q(n+2,1)-q(n+1,1);
    q(n+4,1)=2*q(n+3,1)-q(n+2,1);
    for l=2:n+2
        alphak(l-1,1)=abs(q(l,1)*q(l-1,1))/(abs(q(l-1,1)*q(l,1))+abs(q(l+1,1)*q(l+2,1)));
    end
    Vv(k+1,1)=(1-alphak(k+1,1))*q(k+1,1)+alphak(k+1,1)*q(k+2,1);
    Tv=Vv./abs(Vv);
   
end
end