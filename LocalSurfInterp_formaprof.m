function [U, V, P] = LocalSurfInterp (n, m, Q)

%%LocalSurfInterp
%% Local Surface interpolation
%% through (n+1)x(m+1) points
%%
%%  [U, V, P] = LocalSurfInterp_formaprof (n, m, Q)


%% get ub, r and u direction tangents
total = 0;
ub = zeros (n+1, 1);
td = zeros (n+1, m+1, 3);

for l = 0:m
    %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    Tu0l = ComputeTu0l (l)
    td(0+1, l+1, 0+1) = Tu0l(l+1);
    r(l+1) = 0;
    for k = 1:n
        %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
        Tukl=ComputeTukl(k,l);
        td(k+1, l+1, 0+1) = Tukl(k+1,l+1);
        d = norm (Q(k+1, l+1, :) - Q(k, l+1, :), 2);
        ub(k+1) = ub(k+1) + d;
        r(l+1)  = r(l+1)  + d;
    end
    total = total + r(l+1);
end

for k=1:n
    ub(k+1) = ub(k) + ub(k+1) / total ;
    ub(isnan(ub)) = 0;
end
ub(n+1) = 1;
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
end
vb(m+1) = 1;

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
        a(isnan(a)) = 0 ;
        P(3*i+2,3*k+1)=Q(i+1,k+1)+a*td(i+1,k+1,1);
        P(3*i+3,3*k+1)=Q(i+2,k+1)-a*td(i+2,k+1,1);
    end
end
%Ciclo for per le righe di Bezier control points

for i=0:m
    for k=0:n-1
        b=s(k+1)*delta_vl(k+1)/3;
        b(isnan(b)) = 0 ;
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
        Tu0l(1,l+1)=Vk(1,l+1)/abs(Vk(1,l+1));
        
    end

%% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
    function Tukl = ComputeTukl (k, l)
        Tukl = 0;
        alphak=zeros(n+1,m+1);
        q=zeros(n+1,m+4);
        q1=diff(Q,1,2);
        q(:,3:end-2)=q1;
        q(:,2)=2*q(:,3)-q(:,4);
        q(:,1)=2*q(:,2)-q(:,3);
        q(:,m+3)=2*q(:,m+2)-q(:,m+1);
        q(:,m+4)=2*q(:,m+3)-q(:,m+2);
        alphak(k+1,l+1)=abs(q(k+1,l+2)*q(k+1,l+1))...
            /(abs(q(k+1,l+1)*q(k+1,l+2))+abs(q(k+1,l+3)*q(k+1,l+4)));
        alphak(isnan(alphak))=0;
        Vk(k+1,l+1)=(1-alphak(k+1,l+1))*q(k+1,l+1)+alphak(k+1,l+1)*q(k+1,l+2);
        Tukl(k+1,l+1)=Vk(k+1,l+1)/abs(Vk(k+1,l+1));
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
        Vv(k+1,1)=(1-alphak(k+1))*q(k+1,1)+alphak(k+1,1)*q(k+2,1);
        Tvk0(k+1,1)=Vv(k+1,1)/abs(Vv(k+1,1));
        
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
        Duvkl = zeros(n+1,m+1) ;
        
        %Calcolo Dvkl
        D_u = r.*td(:,:,1);
        
        %Calcolo Dvkl
        
        D_v =s.*td(:,:,2);
        
        %computing delta_uk
        delta_uk = diff(ub') ;
        
        %computing delta_vl
        delta_vl = diff(vb');
        
        
        %Calcolo alpha_k and beta_l
        
        a_k = delta_uk(1:end-1)./(delta_uk(1:end-1)+delta_uk(2:end));
        b_l = delta_vl(1:end-1)./(delta_vl(1:end-1)+delta_vl(2:end));
        q_k = diff(Q,1,2);
        d_k = q_k./delta_uk ;
        
        %Calcolo duvkl e dvukl
        
        d_vu = (1-a_k) .* (D_v(2:end-1,2:end-1)-D_v(1:end-2,2:end-1)) ./(delta_uk(1:end-1)) ...
            + a_k.*(D_v(3:end,2:end-1)-D_v(2:end-1,2:end-1))./ delta_uk(2:end) ;
        
        d_uv= (1-b_l) .* (D_u(2:end-1,2:end-1)-D_u(1:end-2,2:end-1))./delta_vl(1:end-1)+b_l...
            .*(D_u(3:end,2:end-1)-D_u(2:end-1,2:end-1))./delta_vl(2:end);
        
        %Calcolo il "cuore" della matrice tramite la formula (9.59)
        
        Duvkl(2:end-1,2:end-1)=(a_k.*d_uv+b_l.*d_vu)./(a_k+b_l);
        
        % End formulas (9.32) and (9.28) for the boundaries (three point
        % method)
        
        
        %Calcolo elementi interni della prima e ultima riga di D_uv (9.28)
        
        Duvkl(1,2:end-1) = (1-a_k).*d_k(1,1:end-1)+a_k.*d_k(1,2:end);
        Duvkl(n+1,2:end-1) = (1-a_k).*d_k(n+1,1:end-1)+a_k.*d_k(n+1,2:end);
        
        %Calcolo primo e ultimo elemento di ogni riga (9.32)
        Duvkl(:,1) = 2*d_k(:,1)-Duvkl(:,2);
        Duvkl(:,n+1) = 2*d_k(:,n)-Duvkl(:,n+1);
        
        
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
