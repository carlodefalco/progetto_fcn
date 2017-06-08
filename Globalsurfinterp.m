function [U,V,P]=Globalsurfinterp(n,m,Q,p,q)
x=Q(:,:,1);
y=Q(:,:,2);
z=Q(:,:,3);
[uk,vl]=Surfmeshpar1(n,m,Q);
%Compute "U" using eq 9.8
nn=n+p+2;
for j=1:n-p
    summ=0;
    for i=j:j+p-1
            summ=uk(i+1)/p+summ;
    end
     uu(j+p+1)=summ;
end
U=zeros(1,nn);
    for jj=1:nn
        if jj<=p+1
            U(1,jj)=0;
        elseif jj>=nn-p
            U(1,jj)=1;
        else
            U(1,jj)=uu(jj);
        end
    end
%Compute "V" using eq 9.8
mm=m+q+2;
for j=1:m-q
    summ=0;
    for i=j:j+q-1
            summ=vl(i+1)/q+summ;
    end
     vv(j+q+1)=summ;
end

V=zeros(1,mm);
    for jj=1:mm
        if jj<=q+1
            V(1,jj)=0;
        elseif jj>=mm-q
            V(1,jj)=1;
        else
            V(1,jj)=vv(jj);
        end
    end
    
  R=zeros(n+1,m+1,3);
    for l=0:m
       %Do curve interpolation through Q[0][l],...,Q[n][l] 
       [x_R(:,l+1)]=alg1_bis(n,p,x(:,l+1),U,uk,1);
       [y_R(:,l+1)]=alg1_bis(n,p,y(:,l+1),U,uk,1);
       [z_R(:,l+1)]=alg1_bis(n,p,z(:,l+1),U,uk,1);

    end

   for i=0:n
       %Do curve interpolation through Q[0][l],...,Q[n][l] 
       [x_P(i+1,:)]=alg1_bis(m,q,x_R(i+1,:)',V,vl,1);
       [y_P(i+1,:)]=alg1_bis(m,q,y_R(i+1,:)',V,vl,1);
       [z_P(i+1,:)]=alg1_bis(m,q,z_R(i+1,:)',V,vl,1);
   end
P(:,:,1)=x_P;
P(:,:,2)=y_P;
P(:,:,3)=z_P;
end