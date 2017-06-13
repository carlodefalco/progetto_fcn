function [U, V, P,td] = LocalSurfInterp (n, m, Q)

%%LocalSurfInterp
%% Local Surface interpolation
%% through (n+1)x(m+1) points
%%
%%  [U, V, P] = LocalSurfInterp_formaprof (n, m, Q)


%% get ub, r and u direction tangents
total = 0;
ub = zeros (n+1, 1);
td = zeros (n+1, m+1, 3);
Vk = zeros(n+1,m+1);
Vv =zeros(n+1,m+1);
delta_uk = zeros(1,n);
delta_vl = zeros(1,m);
alpha_k = zeros(1,n-1);
beta_l = zeros(1,m-1);
d_vu = zeros(n-1,m-1);

for l = 0:m
    %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    Tu0l = ComputeTu0l (l)
    td(0+1, l+1, 0+1) = Tu0l(l+1);
    r(l+1) = 0;
    for k = 1:n
        %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
        Tukl=ComputeTukl(k,l);
        td(k+1, l+1, 0+1) = Tukl(k,l+1);
        d = norm (Q(k+1, l+1, :) - Q(k, l+1, :), 2);
        ub(k+1) = ub(k+1) + d;
        r(l+1)  = r(l+1)  + d;
    end
    total = total + r(l+1);
end

for k=1:n
    ub(k+1) = ub(k) + ub(k+1) / total ;
    ub(n+1) = 1;
    ub(isnan(ub))=0;
end

%% get vb, s and v direction tangents
total = 0;
vb = zeros (m+1, 1);

for k = 0:n
    %% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)
    Tvk0 = ComputeTvk0 (k) ;
    td(k+1, 0+1, 1+1) = Tvk0(k+1);
    s(k+1) = 0;
    for l = 1:m
        %% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
        Tvkl = ComputeTvkl (k, l)
        td(k+1, l+1, 1+1) = Tvkl(k+1,l+1);
        d = norm (Q(k+1, l+1, :) - Q(k+1, l, :), 2);
        vb(l+1) = vb(l+1) + d;
        s(k+1)  = s(k+1)  + d;
    end
    total = total + s(k+1);
end

for l=1:m
    vb(l+1) = vb(l) + vb(l+1) / total ;
    vb(m+1) = 1;
end

%% TODO: Load the U knot vector
%%       Load the V knot vector
[U, V] = LoadUV ();



%% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k+1, l+1, 2+1)
for k=0:n
    for l=0:m
Duvkl = ComputeDuvkl (k, l);
td(k+1, l+1, 2+1) = Duvkl(k+1,l+1) ;
    end
end


%%       Compute all Bezier control points along each
%%       row and column of data points.

for i=0:n
    for j=0:m
        P(3*j+1,3*i+1)=Q(j+1,i+1);
    end
end

%Ciclo for per le colonne di Bezier control points
for k=0:n
    for i=0:m-1
        
        a=r(k+1)*delta_uk(i+1)/3;
        P(3*i+2,3*k+1)=Q(i+1,k+1)+a*td(i+1,k+1,1);
        P(3*i+3,3*k+1)=Q(i+2,k+1)-a*td(i+2,k+1,1);
    end
end
%Ciclo for per le righe di Bezier control points

for i=0:m
    for k=0:n-1
        b=s(k+1)*delta_vl(k+1)/3;
        P(3*i+1,3*k+2)=Q(i+1,k+1)+b*td(i+1,k+1,2);
        P(3*i+1,3*k+3)=Q(i+1,k+2)-b*td(i+1,k+2,2);
    end
end




%% TODO: Compute the four inner control points of the (k,l)th
%% Bezier patch and load them into P
gamma=(delta_uk*delta_vl.')/9;
for k=0:n-1
    for i=0:m-1
        
        
        P(3*i+2,3*k+2)=gamma*td(i+1,k+1,3)+P(i+1,k+2)+P(i+2,k+1)-P(i+1,k+1);
        P(3*i+3,3*k+2)=-gamma*td(i+2,k+1,3)+P(i+4,k+2)-P(i+4,k+1)+P(i+3,k+1);
        P(3*i+2,3*k+3)=-gamma*td(i+1,k+2,3)+P(i+2,k+4)-P(i+1,k+4)+P(i+1,k+3);
        P(3*i+3,3*k+3)=+gamma*td(i+2,k+2,3)+P(i+3,k+4)-P(i+4,k+3)+P(i+4,k+4);
        
    end
end



%% TODO: Load the NURBS control points by discarding Bezier points
%% along inner rows and columns, Figure 9.32c
P = DiscardInnerRowsandColumns (P);

%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    function Tu0l = ComputeTu0l (l)
        Tu0l = 0;
        alphak=zeros(n+1,m+1);
        q=zeros(n+1,m+4);
        q1=diff(Q,1,2);
        q(:,3:end-2)=q1;
        q(1,2)=2*q(1,3)-q(1,4);
        q(1,1)=2*q(1,2)-q(1,3);
        q(1,m+3)=2*q(1,m+2)-q(1,m+1);
        q(1,m+4)=2*q(1,m+3)-q(1,m+2);
        alphak(1,l+1)=abs(q(1,l+2)*q(1,l+1))/(abs(q(1,l+1)*q(1,l+2))+abs(q(1,l+3)*q(1,l+4)));
        alphak(isnan(alphak))=0;
        Vk(1,l+1)=(1-alphak(1,l+1))*q(1,l+1)+alphak(1,l+1)*q(1,l+2);
        Tu0l=Vk./abs(Vk);
        
    end

%% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
    function Tukl = ComputeTukl (k, l)
        Tukl = 0;
        alphak=zeros(n+1,m+1);
        q=zeros(n+1,m+4);
        q1=diff(Q,1,2);
        q(:,3:end-2)=q1;
        q(k+1,2)=2*q(k+1,3)-q(k+1,4);
        q(k+1,1)=2*q(k+1,2)-q(k+1,3);
        q(k+1,m+3)=2*q(k+1,m+2)-q(k+1,m+1);
        q(k+1,m+4)=2*q(k+1,m+3)-q(k+1,m+2);
        alphak(k+1,l+1)=abs(q(k+1,l+2)*q(k+1,l+1))...
            /(abs(q(k+1,l+1)*q(k+1,l+2))+abs(q(k+1,l+3)*q(k+1,l+4)));
        alphak(isnan(alphak))=0;
        Vk(k+1,l+1)=(1-alphak(k+1,l+1))*q(k+1,l+1)+alphak(k+1,l+1)*q(k+1,l+2);
        Tukl=Vk./abs(Vk); 
    end

%% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)
    function Tvk0 = ComputeTvk0 (k)
        Tvk0 = 0;
        alphak=zeros(n+1,m+1);
        q=zeros(n+4,m+1);
        q1=diff(Q,1,1);
        q(3:end-2,:)=q1;
        q(2,1)=2*q(3,1)-q(4,1);
        q(1,1)=2*q(2,1)-q(3,1);
        q(n+3,1)=2*q(n+2,1)-q(n+1,1);
        q(n+4,1)=2*q(n+3,1)-q(n+2,1);
        alphak(k+1,1)=abs(q(k+2,1)*q(k+1,1))/(abs(q(k+1,1)*q(k+2,1))+abs(q(k+3,1)*q(k+4,1)));
        alphak(isnan(alphak))=0;
        Vv(k+1,1)=(1-alphak(k+1,1))*q(k+1,1)+alphak(k+1,1)*q(k+2,1);
        Tvk0=Vv./abs(Vv);
        
    end

%% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
    function Tvkl = ComputeTvkl (k, l)
        Tvkl = 0;
        alphak=zeros(n+1,m+1);
        q=zeros(n+4,m+1);
        q1=diff(Q,1,1);
        q(3:end-2,:)=q1;
        q(2,l+1)=2*q(3,l+1)-q(4,l+1);
        q(1,l+1)=2*q(2,l+1)-q(3,l+1);
        q(m+3,l+1)=2*q(m+2,l+1)-q(m+1,l+1);
        q(m+4,l+1)=2*q(m+3,l+1)-q(m+2,l+1);
        alphak(k+1,l+1)=abs(q(k+2,l+1)*q(k+1,l+1))/(abs(q(k+1,l+1)*q(k+2,l+1))+abs(q(k+3,l+1)*q(k+4,l+1)));
        alphak(isnan(alphak))=0;
        Vv(k+1,l+1)=(1-alphak(k+1,l+1))*q(k+1,l+1)+alphak(k+1,l+1)*q(k+2,l+1);
        Tvkl=Vv./abs(Vv);
    end


%% TODO: Load the U knot vector
%%       Load the V knot vector
    function [U, V] = LoadUV ()
        U = 0;
        V = 0;
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




%% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k, l, 2)
    function Duvkl = ComputeDuvkl (k, l)
        Duvkl = 0;
        
        for k=0:n
            for l=0:m
                D_u(k+1,l+1)=r(l+1)*td(k+1,l+1,1);
            end
        end
        %Calcolo Dvkl
        for l=0:m
            for k=0:m
                D_v(k+1,l+1)=s(k+1)*td(k+1,l+1,2);
            end
        end
        
        for k=1:n
            delta_uk(k)=ub(k+1)-ub(k);
        end
        %computing delta_vl
        
        for l=1:m
            delta_vl(l)=vb(l+1)-vb(l);
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
        q_k=diff(Q,1,2);
        for k=1:n+1
            for l=1:m
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

%% TODO: Compute the four inner control points of the (k,l)th
%% Bezier patch and load them into P
    function P = ComputeInnerControlPoints ()
        P=zeros(2*n+2,2*m+2,3);
        gamma=(delta_uk*delta_vl.')/9;
        for k=0:n-1
            for i=0:m-1
                
                
                P(3*i+2,3*k+2,:)=gamma*td(i+1,k+1,3)+P(i+1,k+2,:)+P(i+2,k+1,:)-P(i+1,k+1,:);
                P(3*i+3,3*k+2,:)=-gamma*td(i+2,k+1,3)+P(i+4,k+2,:)-P(i+4,k+1,:)+P(i+3,k+1,:);
                P(3*i+2,3*k+3,:)=-gamma*td(i+1,k+2,3)+P(i+2,k+4,:)-P(i+1,k+4,:)+P(i+1,k+3,:);
                P(3*i+3,3*k+3,:)=+gamma*td(i+2,k+2,3)+P(i+3,k+4,:)-P(i+4,k+3,:)+P(i+4,k+4,:);
                
            end
        end
    end

%% TODO: Load the NURBS control points by discarding Bezier points
%% along inner rows and columns, Figure 9.32c
    function P = DiscardInnerRowsandColumns (P)
        for i=m-2:-1:0
            P(3*i+4,:,:)=[];
        end
        
        for i=m-2:-1:0
            P(:,3*i+4,:)=[];
        end
    end

end
