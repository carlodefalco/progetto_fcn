function [U,V]=LoadUV(ub,vb,n,m)
%Load the U knot vector (9.56)
e=8+2*(n-1);
w=8+2*(m-1);
U=zeros(1,e);
V=zeros(1,e);

for j=1:4
    U(1,j)=0;
end
for j=e-4:e
    U(1,j)=1;
end
for j=4:e-4
    a=floor((j-1)/2);
    U(1,j)=ub(a);
end

%Load the V knot vector (9.56)
for j=1:4
    V(1,j)=0;
end
for j=w-4:w
    V(1,j)=1;
end
for j=4:w-4
    a=floor((j-1)/2);
    V(1,j)=vb(a);
end
end