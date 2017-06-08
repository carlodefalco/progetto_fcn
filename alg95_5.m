clc
clear all
out = A93_94_95_test_data (100 , 100, 17, 25);
Q(:,:,1)=out{1,1,1};
Q(:,:,2)=out{1,1,2};
Q(:,:,3)=out{1,1,3};

[n, m, o]=size(Q);
n=n-1;
m=m-1;
x=Q(:,:,1);
y=Q(:,:,2);
z=Q(:,:,3);




total=0;
ub=zeros(1,n+1);
q=zeros(n+1,n+4);
r=zeros(1,n+1);
%get ub,r and direction tangents
for k=0:n
    ub(k+1)=0;
end

for l=0:m
    
    %Ciclo for per calcolare la matrice delle distanze tra i punti
    for k=0:n
        for i=0:m-1
            %Compute and load T0,l,u into td_0[0][l]
            q(k+1,i+3)=sqrt((x(k+1,i+2)-x(k+1,i+1))^2+(y(k+1,i+2)-y(k+1,i+1))^2+(z(k+1,i+2)-z(k+1,i+1))^2);
            
        end
        
        q(k+1,2)=2*q(k+1,3)-q(k+1,4);
        q(k+1,1)=2*q(k+1,2)-q(k+1,3);
        q(k+1,n+3)=2*q(k+1,n+2)-q(k+1,n+1);
        q(k+1,n+4)=2*q(k+1,n+3)-q(k+1,n+2);
    end
    r(l+1)=0;
    td=zeros(n+1,m+1,3);
    %Ciclo for per calcolare alpha_k
    %     %Ciclo for per calcolare V_k e caricare T_k in td_0[0][l]
    for k=2:n+1
        alpha(1,k)=abs(q(1,k)*q(1,k-1))/(abs(q(1,k-1)*q(1,k))+abs(q(1,k+1)*q(1,k+2)));
    end
    for k=1:n+1
        Vu(1,k)=(1-alpha(1,k))*q(1,k)+alpha(1,k)*q(1,k+1);
        T_ku(1,k)=Vu(1,k)/abs(Vu(1,k));
    end
    
    
    for k=1:n
        %Compute and load Tk,l,u into td_0[k][l]
        %Ciclo for per calcolare alpha_k
        for ll=2:n+1
            alpha(k+1,ll)=abs(q(k+1,ll-1)*q(k+1,ll))/(abs(q(k+1,ll-1)*q(k+1,ll))+abs(q(k+1,ll+1)*q(k+1,ll+2)));
            %         end
            %         %Ciclo for per calcolare V_k e caricare T_k in td_0[k][l]
            %         for l=1:m+1
        end
        for ll=1:n+1
            Vu(k+1,ll)=(1-alpha(k+1,ll))*q(k+1,ll)+alpha(k+1,ll)*q(k+1,ll+1);
            T_ku(k+1,ll)=Vu(k+1,ll)/abs(Vu(k+1,ll));
        end
        
        
        d=sqrt((x(k+1,l+1)-x(k,l+1))^2+(y(k+1,l+1)-y(k,l+1))^2+(z(k+1,l+1)-z(k,l+1))^2);
        ub(k+1)=ub(k+1)+d;
        r(l+1)=r(l+1)+d;
        
    end
    
    total=total+r(l+1);
end
td(:,:,1)=T_ku;



%Ciclo for per calcolare ub
for k=1:n
    ub(k+1)=ub(k)+ub(k+1)/total;
end

%get vb,s and direction tangents
total=0;
q=zeros(n+1,m+4);
s=zeros(1,n+1);
vb=zeros(1,n+1);
alpha=zeros(1,n+1);
for l=0:m
    vb(l+1)=0;
end
for k=0:n
    %Ciclo for per calcolare la matrice delle distanze tra i punti
    for l=0:m
        for i=0:n-1
            %Compute and load T0,l,u into td_0[0][l]
            q(l+1,i+3)=sqrt((x(i+2,l+1)-x(i+1,l+1))^2+(y(i+2,l+1)-y(i+1,l+1))^2+(z(i+2,l+1)-z(i+1,l+1))^2);
            
        end
        
        q(l+1,2)=2*q(l+1,3)-q(l+1,4);
        q(l+1,1)=2*q(l+1,2)-q(l+1,3);
        q(l+1,n+3)=2*q(l+1,n+2)-q(l+1,n+1);
        q(l+1,n+4)=2*q(l+1,n+3)-q(l+1,n+2);
    end
    s(k+1)=0;
    %Ciclo for per calcolare alpha_k
    %Since q is already changed wrt the previous case there's no need for
    %changing alpha
    for l=2:m+1
        alpha(l,1)=abs(q(1,l-1)*q(1,l))/(abs(q(1,l-1)*q(1,l))+abs(q(1,l+1)*q(1,l+2)));
    end
    for l=1:m+1
        Vv(l,1)=(1-alpha(l,1))*q(1,l)+alpha(l,1)*q(1,l+1);
        T_kv(l,1)=Vv(l,1)/abs(Vv(l,1));
    end
    
    
    
    for l=1:m
        %Compute and load Tk,l,u into td_1[k][l]
        %Ciclo for per calcolare alpha_k
        for kk=2:m+1
            alpha(kk,l+1)=abs(q(l+1,kk-1)*q(l+1,kk))/(abs(q(l+1,kk-1)*q(l+1,kk))+abs(q(l+1,kk+1)*q(l+1,kk+2)));
        end
        %         %Ciclo for per calcolare V_k e caricare T_k in td_0[k][l]
        for kk=1:m+1
            Vv(kk,l+1)=(1-alpha(kk,l+1))*q(l+1,kk)+alpha(kk,l+1)*q(l+1,kk+1);
            T_kv(kk,l+1)=Vv(kk,l+1)/abs(Vv(kk,l+1));
        end
        
        
        d1=sqrt((x(k+1,l+1)-x(k+1,l))^2+(y(k+1,l+1)-y(k+1,l))^2+(z(k+1,l+1)-z(k+1,l))^2);
        
        vb(l+1)=vb(l+1)+d1;
        s(k+1)=s(k+1)+d1;
    end
    total=total+s(k+1);
    
end
td(:,:,2)=T_kv;

%Ciclo for per calcolare vb
for l=1:m
    vb(l+1)=vb(l)+vb(l+1)/total;
    vb(m+1)=1;
end



%Load the U knot vector (9.56)
e=8+2*(n-1);
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
    U(1,j)=ub(1,a);
end

%Load the V knot vector (9.56)
for j=1:4
    V(1,j)=0;
end
for j=e-4:e
    V(1,j)=1;
end
for j=4:e-4
    a=floor((j-1)/2);
    V(1,j)=vb(1,a);
end








% delta_uk=zeros(1,n);
%computing delta_uk
for k=1:n
    delta_uk(k)=ub(k+1)-ub(k);
end
%computing delta_vl

for l=1:m
    delta_vl(l)=vb(l+1)-vb(l);
end
%Compute all Beziér control points along each row and column of data
%points
%Calcolo matrice punti di controllo (andrà messa in mezzo al codice)

%Ciclo for per inserire control points estremi
for i=0:n
    for j=0:m
        P(3*j+1,3*i+1,:)=Q(j+1,i+1,:);
    end
end

%Ciclo for per le colonne di Bezier control points
for k=0:n
    for i=0:m-1
        
        a=r(i+1)*delta_uk(i+1)/3;
        P(3*i+2,3*k+1,:)=Q(i+1,k+1,:)+a*T_ku(i+1,k+1);
        P(3*i+3,3*k+1,:)=Q(i+2,k+1,:)-a*T_ku(i+2,k+1);
    end
end
%Ciclo for per le righe di Bezier control points

for i=0:m
    for k=0:n-1
        b=s(k+1)*delta_vl(k+1)/3;
        P(3*i+1,3*k+2,:)=Q(i+1,k+1,:)+b*T_kv(i+1,k+1);
        P(3*i+1,3*k+3,:)=Q(i+1,k+2,:)-b*T_kv(i+1,k+2);
    end
end



%Calcolo D_u e D_v

for k=1:n+1
    for l=1:m+1
        D_u(k,l)=r(l)*T_ku(k,l);
    end
end

for l=1:m+1
    for k=1:m+1
        D_v(k,l)=s(k)*T_kv(k,l);
    end
end



%computing alpha_k and beta_l

for k=1:n-1
    alpha_k(k)=delta_uk(k)/(delta_uk(k)+delta_uk(k+1));
end


for l=1:m-1
    beta_l(l)=delta_vl(l)/(delta_vl(l)+delta_vl(l+1));
end

%Calcolo d_uv e d_vu

for l=1:m-1
    for k=1:n-1
        d_vu(k,l)=(1-alpha_k(k))*(D_v(k+1,l)-D_v(k,l))/delta_uk(k)+alpha_k(1,k)*(D_v(k+2,l)-D_v(k+1,l))/delta_uk(k+1) ;
        
    end
end
for k=1:n-1
    for l=1:m-1
        d_uv(k,l)=(1-beta_l(1,l))*(D_u(k,l+1)-D_u(k,l))/delta_vl(1,l)+beta_l(1,l)*(D_u(k,l+2)-D_u(k,l+1))/delta_vl(1,l+1);
    end
end

%Calcolo il "cuore" della matrice D_uv, i boundaries verranno calcolati
%successivamente
for k=1:n-1
    for l=1:m-1
        D_uv(k+1,l+1)=(alpha_k(k)*d_uv(k,l)+beta_l(l)*d_vu(k,l))/(alpha_k(k)+beta_l(l));
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
    D_uv(1,k)=(1-alpha_k(k-1))*d_k(1,k-1)+alpha_k(k-1)*d_k(1,k);
    D_uv(n+1,k)=(1-alpha_k(k-1))*d_k(n+1,k-1)+alpha_k(k-1)*d_k(n+1,k);
end

%Calcolo primo e ultimo elemento di ogni riga
for l=1:n+1
    D_uv(l,1)=2*d_k(l,1)-D_uv(l,2);
    D_uv(l,n+1)=2*d_k(l,n)-D_uv(l,n-1);
end
td(:,:,3)=D_uv;


%Compute the four inner control points of the (k,l)th Bezier patch and load
%them into P

gamma=(delta_uk*delta_vl.')/9;

for k=0:n-1
    for i=0:m-1
        
        
        P(3*i+2,3*k+2,:)=gamma*D_uv(i+1,k+1)+P(i+1,k+2,:)+P(i+2,k+1,:)-P(i+1,k+1,:);
        P(3*i+3,3*k+2,:)=-gamma*D_uv(i+2,k+1)+P(i+4,k+2,:)-P(i+4,k+1,:)+P(i+3,k+1,:);
        P(3*i+2,3*k+3,:)=-gamma*D_uv(i+1,k+2)+P(i+2,k+4,:)-P(i+1,k+4,:)+P(i+1,k+3,:);
        P(3*i+3,3*k+3,:)=+gamma*D_uv(i+2,k+2)+P(i+3,k+4,:)-P(i+4,k+3,:)+P(i+4,k+4,:);
        
    end
end


%Load the Nurbs control points by discarding Bézier points along inner rows
%and columns
for i=m-2:-1:0
    P(3*i+4,:,:)=[];
end

for i=m-2:-1:0
    P(:,3*i+4,:)=[];
end



knots = {U V} ;
% cntrl = zeros [3 , n , m];
cntrl(1,:,:)=P(:,:,1);
cntrl(2,:,:)=P(:,:,2);
cntrl(3,:,:)=P(:,:,3);
nrb = nrbmak(cntrl,knots);
nrbplot(nrb,[100 100]);






