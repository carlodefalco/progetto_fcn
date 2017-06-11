function [Tu]=ComputeTukl(n,m,Q,l,k)
for l=0:m
    for k=1:n
        alphak=zeros(n+1,m+1);
        q=zeros(n+1,m+4);
        for i=1:m
            q(k+1,i+2)=sqrt((Q(k+1,i+1,1)-Q(k+1,i,1))^2+(Q(k+1,i+1,2)-Q(k+1,i,2))^2+(Q(k+1,i+1,3)-Q(k+1,i,3))^2);
        end
        q(k+1,2)=2*q(k+1,3)-q(k+1,4);
        q(k+1,1)=2*q(k+1,2)-q(k+1,3);
        q(k+1,m+3)=2*q(k+1,m+2)-q(k+1,m+1);
        q(k+1,m+4)=2*q(k+1,m+3)-q(k+1,m+2);
        for i=2:m+2
            alphak(k+1,i-1)=abs(q(k+1,i)*q(k+1,i-1))/(abs(q(k+1,i-1)*q(k+1,i))+abs(q(k+1,i+1)*q(k+1,i+2)));
        end
        Vk(k+1,l+1)=(1-alphak(k+1,l+1))*q(k+1,l+1)+alphak(k+1,l+1)*q(k+1,l+2);
        Tu=Vk./abs(Vk);
    end
end
end