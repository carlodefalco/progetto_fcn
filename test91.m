function [U,m,P]=alg91(n,p,Q,r)
%Q=[x;y]';
x=Q(:,1);
y=Q(:,2);
z=Q(:,3);
n=numel(x)-1;
%Calcolo distanze tra punti
for i=1:n
d(i)=abs([x(i+1)-x(i)]^2+[y(i+1)-y(i)]^2+[z(i+1)-z(i)]^2)^(1/2);
end
d_tot=sum(d);

%Calcolo dei parametri u segnato
u=zeros(n+1,1)';
u(1)=0;
u(numel(x))=1;

for i=1:(n-1)
u(i+1)=u(i)+d(i)/d_tot;
end

%Calcolo knot vector "U" usando eq 9.8
m=n+p+2; 
for j=1:n-p
    summ=0;
    for i=j:j+p-1
            summ=u(i+1)/p+summ;
    end
     uu(j+p+1)=summ;
end

U=zeros(1,m);
    for jj=1:m
        if jj<=p+1
            U(1,jj)=0;
        elseif jj>=m-p
            U(1,jj)=1;
        else
            U(1,jj)=uu(jj);
        end
    end
    
%Initialize array A to zero
A=zeros(n+1);
%Set up coefficient matrix
% Bisogna aggiungere una condizione per calcolare la matrice, ossia che
% N(i,p)(u(k))=0 se |i-k|>=p, in questo modo dovrebbe uscire la matrice di
% dimensione giusta.
%Vedere propriet? 2.1
for ii=1:n+1
    N=zeros(1,n+1);
    span=findspan(n+1,p,u(ii),U);
    
for j=1:n+1
        if span<j  ||  span>=j+p+1
           N(1,j)=0;
        else          
            N(1,span-p+1:span+1)=basisfuncs(span,u(ii),p,U);
            A(ii,:)=N;
        end 
end 
end
A

% [L,UU]=lu(A);
% xx_P=Q(:,1);
% yy_P=Q(:,2);
% x_P=avanti(L,xx_P);
% y_P=avanti(L,yy_P);
% x=indietro(UU,x_P);
% y=indietro(UU,y_P);
% % y=avanti(L,Q);
% % PP=avanti(U,y);
% %PP=A\Q;
% P=[x y];

[L,UU,PP]=lu(A);
y=L\(PP*Q);
P=UU\y;


 
%verifica
%P=A\Q;
