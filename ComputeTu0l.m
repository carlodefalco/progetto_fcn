function [Tu]=ComputeTu0l(n,m,Q,l)
for l=0:m
    alphak=zeros(n+1,m+1);
    q=zeros(n+1,m+4);
    for i=1:m
        q(1,i+2)=sqrt((Q(1,i+1,1)-Q(1,i,1))^2+(Q(1,i+1,2)-Q(1,i,2))^2+(Q(1,i+1,3)-Q(1,i,3))^2);
    end
    q(1,2)=2*q(1,3)-q(1,4);
    q(1,1)=2*q(1,2)-q(1,3);
    q(1,m+3)=2*q(1,m+2)-q(1,m+1);
    q(1,m+4)=2*q(1,m+3)-q(1,m+2);
    for k=2:m+2
        alphak(1,k-1)=abs(q(1,k)*q(1,k-1))/(abs(q(1,k-1)*q(1,k))+abs(q(1,k+1)*q(1,k+2)));
    end
    Vk(1,l+1)=(1-alphak(1,l+1))*q(1,l+1)+alphak(1,l+1)*q(1,l+2);
    Tu=Vk./abs(Vk);
   
end
end