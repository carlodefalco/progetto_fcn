function [U,P]=alg96(Q,r,Wq,D,s,Wd,n,p)
Q=zeros(n,3);
x=Q(:,1);
y=Q(:,2);
z=Q(:,3);
ru=-1;
rc=-1;
for i=0:r
    if Wq(i+1)>0
ru=ru+1;
    else
        rc=rc+1;
    end
end
su=-1;
sc=-1;
for j=0:s
    if Wd(j+1)>0
        su=su+1;
    else
        sc=sc+1;
    end
end
mu=ru+su+1;
mc=rc+sc+1;
if mc>=n || mc+n>=mu+1
    display('Error')
end

%Compute and load parameters u_k into ub[] (Eq.[9.5]);
for i=1:r
dd(i)=abs([x(i+1)-x(i)]^2+[y(i+1)-y(i)]^2+[z(i+1)-z(i)]^2)^(1/2);
end
d_tot=sum(dd);
u_k=zeros(r+1,1)';
u_k(1)=0;
u_k(numel(x))=1;
for i=1:(r-1)
u_k(i+1)=u_k(i)+dd(i)/d_tot;
end
ub=u_k;

%Compute and load the knots into U[] (Eqs.[9.68],[9.69]);
d=(r+1)/(n-p+1);
nn=n+p+2;
U=zeros(1,nn);

for j=1:n-p
    i=ceil(j*d);
    alpha=j*d-i;
    U(p+j+1)=(1-alpha)*ub(i)+alpha*ub(i+1);
end
% for j=1:n-p
%     summ=0;
%     for i=j:j+p-1
%             summ=ub(i+1)/p+summ;
%     end
%      uu(j+p+1)=summ;
% end

    for jj=1:nn
        if jj<=p+1
            U(1,jj)=0;
        elseif jj>=nn-p
            U(1,jj)=1;
%         else
%             U(1,jj)=uu(jj);
        end
    end
   %Now set up arrays N,W,S,T,M
   j=0;
   mu2=0;
   mc2=0;
   N=zeros(r+1,n+1);
   for i=0:r
       span=findspan(n+1,p,ub(i+1),U);
       ss(i+1)=span;
       dflag=0;
       if j<=s
           if i==I(j)
               dflag=1;
           end
       end
           if dflag==0
               
            funs=zeros(1,n+1);
            funs(1,span-p+1:span+1)=basisfuncs(span,ub(i+1),p,U);
           else
               funs=DersBasisFuns(span,ub(i+1),p,1,U);
           end
           if Wq(i+1)>0                %Unconstrained points

               W(mu2+1)=Wq(i+1);
               %Load the mu2th row of N[][] from funs[0][];
               N(mu2+1,:)=funs;
               x_S(mu2+1)=W(mu2+1)*x(i+1);
               y_S(mu2+1)=W(mu2+1)*y(i+1);
               z_S(mu2+1)=W(mu2+1)*z(i+1);
               mu2=mu2+1;
           else                        %Constrained points
   
               %Load the mc2th row of M[][] from funs [0][];
               M(mc2+1,:)=funs;
               T(mc2+1)=Q(i+1);
               mc2=mc2+1;
             end  
   end
  W=diag(W);
  
   NN=N'*W*N;
   x_SS=N'*x_S';
   y_SS=N'*y_S';
   z_SS=N'*z_S';
   [L,UU]=lu(NN);
   xx_SS=avanti(L,x_SS);
   yy_SS=avanti(L,y_SS);
   zz_SS=avanti(L,z_SS);
   %NN_inv=inv(NN);
   if mc<0
       x_P=avanti(UU,xx_SS);
       y_P=avanti(UU,yy_SS);
       z_P=avanti(UU,zz_SS);
   end
   
   P=[x_P y_P z_P];
end
 
