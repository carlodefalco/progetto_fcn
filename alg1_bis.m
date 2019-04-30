function [P]=alg1_bis(n,p,Q,U,u,r)

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

[L,UU]=lu(A);
xx_P=Q(:,1);
x_P=avanti(L,xx_P);
x=indietro(UU,x_P);
P=x;

