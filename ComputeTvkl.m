function [Tv]=ComputeTvkl(n,m,Q,k,l)
for k=0:n
    for l=1:m
        alphak=zeros(n+1,m+1);
        q=zeros(n+4,m+1);
        for i=1:n
            q(i+2,l+1)=sqrt((Q(i+1,l+1,1)-Q(i+1,l+1,1))^2+(Q(i+1,l+1,2)-Q(i,l+1,2))^2+(Q(i+1,l+1,3)-Q(i,l+1,3))^2);
        end
        q(2,l+1)=2*q(3,l+1)-q(4,l+1);
        q(1,l+1)=2*q(2,l+1)-q(3,l+1);
        q(m+3,l+1)=2*q(m+2,l+1)-q(m+1,l+1);
        q(m+4,l+1)=2*q(m+3,l+1)-q(m+2,l+1);
        for i=2:n+2
            alphak(i-1,l+1)=abs(q(i,l+1)*q(i-1,l+1))/(abs(q(i-1,l+1)*q(i,l+1))+abs(q(i+1,l+1)*q(i+2,l+1)));
        end
        Vv(k+1,l+1)=(1-alphak(k+1,l+1))*q(k+1,l+1)+alphak(k+1,l+1)*q(k+2,l+1);
        Tv=Vv./abs(Vv);
    end
end
end