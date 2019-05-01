function x=avanti(A,b)
n=length(b);
x=zeros(n,1);
x(1)=b(1)/A(1,1);
for k=2:n
    x(k)=b(k)-A(k,1:k-1)*x(1:k-1);
    x(k)=x(k)/A(k,k);
end
end
