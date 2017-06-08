function [U,V,P]=GlobalSurfApproximation(r,s,Q,p,q,n,m)
%Global surface approx with fuxed number of control points

x=Q(:,:,1);
y=Q(:,:,2);
z=Q(:,:,3);
[ub,vb]=Surfmeshpar1(r,s,Q);


% Compute knots U by Eqs. (9.68),(9.69);
d=(r+1)/(n-p+1);
nn=n+p+2;
U=zeros(1,nn);

for j=1:n-p
    i=ceil(j*d);
    alpha=j*d-i;
    U(p+j+1)=(1-alpha)*ub(i)+alpha*ub(i+1);
end
for jj=1:nn
        if jj<=p+1
            U(1,jj)=0;
        elseif jj>=nn-p
            U(1,jj)=1;
        end
end
% Compute knots V by Eqs. (9.68),(9.69);
d=(s+1)/(m-q+1);
mm=m+q+2;
V=zeros(1,mm);

for j=1:m-q
    i=ceil(j*d);
    alpha=j*d-i;
    V(p+j+1)=(1-alpha)*vb(i)+alpha*vb(i+1);
end
for jj=1:mm
        if jj<=q+1
            V(1,jj)=0;
        elseif jj>=mm-q
            V(1,jj)=1;
        end
end

%Compute Nu[][] and NTNu[][] using Eq.(9.66);
   for i=0:r
       span=findspan(n,p,ub(i+1),U);
       funs=zeros(1,n+1);
            funs(1,span-p+1:span+1)=basisfuncs(span,ub(i+1),p,U);
         N(i+1,:)=funs;    
   end
  
   for jj=1:r-1
   Nu(jj,:)=N(jj+1,2:n);
   end
      NTNu=Nu'*Nu;

%LUDecomposition(NTNu,n-1,p);
[L,UU,M]=lu(NTNu);
Temp=zeros(n+1,s+1);
for j=0:s
    Temp(1,j+1,1)=Q(1,j+1,1);
    Temp(1,j+1,2)=Q(1,j+1,2);
    Temp(1,j+1,3)=Q(1,j+1,3);

    Temp(n+1,j+1,1)=Q(r+1,j+1,1);
    Temp(n+1,j+1,2)=Q(r+1,j+1,2);
    Temp(n+1,j+1,3)=Q(r+1,j+1,3);

    %Compute and load Ru[] (Eqs.[9.63],[9.67]);

    for k=1:r-1
        RRu(k,j+1,1)=Q(k+1,j+1,1)-N(k+1,1)*Q(1,j+1,1)-N(k+1,n+1)*Q(r+1,j+1,1);
        RRu(k,j+1,2)=Q(k+1,j+1,2)-N(k+1,1)*Q(1,j+1,2)-N(k+1,n+1)*Q(r+1,j+1,2);
        RRu(k,j+1,3)=Q(k+1,j+1,3)-N(k+1,1)*Q(1,j+1,3)-N(k+1,n+1)*Q(r+1,j+1,3);
    end
    
    %Here you get matrix as in Eq. 9.67
    Ru(:,j+1,1)=Nu'*RRu(:,j+1,1);
    Ru(:,j+1,2)=Nu'*RRu(:,j+1,2);
    Ru(:,j+1,3)=Nu'*RRu(:,j+1,3);
    
    
    Temp1(:,j+1,1)=NTNu\Ru(:,j+1,1);
    Temp1(:,j+1,2)=NTNu\Ru(:,j+1,2);
    Temp1(:,j+1,3)=NTNu\Ru(:,j+1,3);

end

for k=1:n-1
    Temp(k+1,:,1)=Temp1(k,:,1);
    Temp(k+1,:,2)=Temp1(k,:,2);
    Temp(k+1,:,3)=Temp1(k,:,3);
end

%Compute Nv[][] and NTNv[][] using Eq.(9.66);
for i=0:s
       span=findspan(m,q,vb(i+1),V);
       funs=zeros(1,m+1);
            funs(1,span-q+1:span+1)=basisfuncs(span,vb(i+1),q,V);
         NN(i+1,:)=funs;    
   end
  
   for jj=1:s-1
   Nv(jj,:)=NN(jj+1,2:m);
   end
      NTNv=Nv'*Nv;
      
%LUDecomposition(NTNv,m-1,q);
[L,UU,M]=lu(NTNv);

for i=0:n
    % v direction fits
    P(i+1,1,1)=Temp(i+1,1,1);
    P(i+1,1,2)=Temp(i+1,1,2);
    P(i+1,1,3)=Temp(i+1,1,3);
    
    P(i+1,m+1,1)=Temp(i+1,s+1,1);
    P(i+1,m+1,2)=Temp(i+1,s+1,2);
    P(i+1,m+1,3)=Temp(i+1,s+1,3);
    %Compute and load Rv[] (Eq.[9.63],[9.67]);
    for k=1:s-1
        RRv(i+1,k,1)=Q(i+1,k+1,1)-NN(k+1,1)*Q(i+1,1,1)-NN(k+1,m+1)*Q(i+1,s+1,1);
        RRv(i+1,k,2)=Q(i+1,k+1,2)-NN(k+1,1)*Q(i+1,1,2)-NN(k+1,m+1)*Q(i+1,s+1,2);
        RRv(i+1,k,3)=Q(i+1,k+1,3)-NN(k+1,1)*Q(i+1,1,3)-NN(k+1,m+1)*Q(i+1,s+1,3);
    end
    
    Rv(i+1,:,1)=Nv'*RRv(i+1,:,1)';
    Rv(i+1,:,2)=Nv'*RRv(i+1,:,2)';
    Rv(i+1,:,3)=Nv'*RRv(i+1,:,3)';
    %Call FrowardBackward() to get the control points
     
         P1(i+1,:,1)=NTNv'\Rv(i+1,:,1)';
         P1(i+1,:,2)=NTNv'\Rv(i+1,:,2)';
         P1(i+1,:,3)=NTNv'\Rv(i+1,:,3)';
end

for k=1:m-1
    P(:,k+1,1)=P1(:,k,1);
    P(:,k+1,2)=P1(:,k,2);
    P(:,k+1,3)=P1(:,k,3);
end
end