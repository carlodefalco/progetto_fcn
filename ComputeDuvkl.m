function [Duvkl]=ComputeDuvkl(n,m,Q,r,s,td,delta_uk,delta_vl)
%Calcolo Dukl

for k=1:n+1
    for l=1:m+1
        D_u(k,l)=r(l)*td(k,l,1);
    end
end
%Calcolo Dvkl
for l=1:m+1
    for k=1:m+1
        D_v(k,l)=s(k)*td(k,l,2);
    end
end



%Calcolo alpha_k and beta_l

for k=1:n-1
    alpha_k(k)=delta_uk(k)/(delta_uk(k)+delta_uk(k+1));
end


for l=1:m-1
    beta_l(l)=delta_vl(l)/(delta_vl(l)+delta_vl(l+1));
end

%Calcolo duvkl e dvukl

for l=1:m-1
    for k=1:n-1
        d_vu(k,l)=(1-alpha_k(k))*(D_v(k+1,l)-D_v(k,l))/delta_uk(k)+alpha_k(k)*(D_v(k+2,l)-D_v(k+1,l))/delta_uk(k+1) ;
        
    end
end
for k=1:n-1
    for l=1:m-1
        d_uv(k,l)=(1-beta_l(l))*(D_u(k,l+1)-D_u(k,l))/delta_vl(l)+beta_l(l)*(D_u(k,l+2)-D_u(k,l+1))/delta_vl(l+1);
    end
end

%Calcolo il "cuore" della matrice D_uv, i boundaries verranno calcolati
%successivamente
for k=1:n-1
    for l=1:m-1
        Duvkl(k+1,l+1)=(alpha_k(k)*d_uv(k,l)+beta_l(l)*d_vu(k,l))/(alpha_k(k)+beta_l(l));
    end
end

%Ciclo for per calcolare gli estremi di D_uv tramite il three point scheme
%(9.32) (9.28)
%Calcolo prima d_k e q_k

for k=1:n+1
    for l=1:m
        q_k(k,l)=sqrt((Q(k,l+1,1)-Q(k,l,1))^2+(Q(k,l+1,2)-Q(k,l,2))^2+(Q(k,l+1,3)-Q(k,l,3))^2);
        d_k(k,l)=q_k(k,l)/delta_uk(l);
    end
end

%Calcolo elementi interni della prima e ultima riga di D_uv
for k=2:n
    Duvkl(1,k)=(1-alpha_k(k-1))*d_k(1,k-1)+alpha_k(k-1)*d_k(1,k);
    Duvkl(n+1,k)=(1-alpha_k(k-1))*d_k(n+1,k-1)+alpha_k(k-1)*d_k(n+1,k);
end

%Calcolo primo e ultimo elemento di ogni riga
for l=1:n+1
    Duvkl(l,1)=2*d_k(l,1)-Duvkl(l,2);
    Duvkl(l,n+1)=2*d_k(l,n)-Duvkl(l,n-1);
end
end