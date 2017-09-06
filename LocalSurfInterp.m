function [U, V, P] = LocalSurfInterp (n, m, Q)

%%LocalSurfInterp
%% Local Surface interpolation
%% through (n+1)x(m+1) points
%%
%%  [U, V, P] = LocalSurfInterp (n, m, Q)


%% get ub, r and u direction tangents
total = 0;
ub = zeros (n+1, 1);
td = zeros (n+1, m+1, 3, 3);

for l = 0:m
    %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    td(0+1, l+1, :, 0+1) = ComputeTu0l (l);
    r(l+1) = 0;
    for k = 1:n
        %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
        td(k+1, l+1, :, 0+1) = ComputeTukl(k,l);
        d = sqrt( (Q(k+1, l+1, 1) - Q(k, l+1, 1))^2+...
            (Q(k+1, l+1, 2) - Q(k, l+1, 2))^2+(Q(k+1, l+1, 3) - Q(k, l+1, 3))^2);
        ub(k+1) = ub(k+1) + d;
        r(l+1)  = r(l+1)  + d;
    end
    total = total + r(l+1);
end

for k=1:n-1
    ub(k+1) = ub(k) + ub(k+1) / total ;
    ub(isnan(ub)) = 0;
end
ub(n+1) = 1;

%% get vb, s and v direction tangents
total = 0;
vb = zeros (m+1, 1);

for k = 0:n
    %% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)
    td(k+1, 0+1, :, 1+1) = ComputeTvk0 (k) ;
    s(k+1) = 0;
    for l = 1:m
        %% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
        td(k+1, l+1, :, 1+1) = ComputeTvkl (k, l) ;
        d = sqrt( (Q(k+1, l+1, 1) - Q(k+1, l, 1))^2+(Q(k+1, l+1, 2) - Q(k+1, l, 2))^2+ ...
            (Q(k+1, l+1, 3) - Q(k+1, l, 3))^2);
        vb(l+1) = vb(l+1) + d;
        s(k+1)  = s(k+1)  + d;
    end
    total = total + s(k+1);
end

for l=1:m-1
    vb(l+1) = vb(l) + vb(l+1) / total ;
end
vb(m+1) = 1;

%% TODO: Load the U knot vector
%%       Load the V knot vector
[U, V] = LoadUV ();


%%       Compute all Bezier control points along each
%%       row and column of data points.
P = zeros(n+1,m+1,3);
for k=0:n
    for l=0:m
        P(3*l+1,3*k+1,:) = Q(l+1,k+1,:);
    end
end

%Ciclo for per le colonne di Bezier control points
for l=0:m
    for k=0:n-1
        delta_uk = diff(ub) ;
        a = r(l+1) * delta_uk(k+1)/3;
        a(isnan(a)) = 0 ;
        td(isnan(td)) = 0 ;
        P(3*k+2,3*l+1,:) = Q(k+1,l+1,:)+a*td(k+1,l+1,:,1);
        P(3*k+3,3*l+1,:) = Q(k+2,l+1,:)-a*td(k+2,l+1,:,1);
    end
end
%Ciclo for per le righe di Bezier control points

for k=0:n
    for l=0:m-1
        delta_vl = diff(vb) ;
        b = s(k+1) * delta_vl(l+1)/3;
        b(isnan(b)) = 0 ;
        P(3*k+1,3*l+2,:) = Q(k+1,l+1,:)+b*td(k+1,l+1,:,2);
        P(3*k+1,3*l+3,:) = Q(k+1,l+2,:)-b*td(k+1,l+2,:,2);
    end
end

%% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k+1, l+1, 2+1)

for k=0:n
    for l=0:m
        td(k+1, l+1, :, 2+1) = ComputeDuvkl (k, l);
    end
end




%% TODO: Compute the four inner control points of the (k,l)th
%% Bezier patch and load them into P
P = ComputeInnerControlPoints (P);


%% TODO: Load the NURBS control points by discarding Bezier points
%% along inner rows and columns, Figure 9.32c
P = DiscardInnerRowsandColumns (P);

%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    function Tu0l = ComputeTu0l (l)
        if ( l < 2 )
            q(1,3:4,:) = diff(Q(1,l+1:l+3,:),1,2);
            q(1,2,:) = 2 * q(1,3,:) - q(1,4,:);
            q(1,1,:) = 2 * q(1,2,:) - q(1,3,:);
        elseif ( l> m-2 )
            q(1,1:2,:) = diff(Q(1,l-2:l,:),1,2);
            q(1,3,:) = 2 * q(1,2,:) - q(1,1,:);
            q(1,4,:) = 2 * q(1,3,:) - q(1,2,:);
        else
            q = diff(Q(1,l-1:l+3,:),1,2);
        end
        c1 = cross(q(1,1,:),q(1,2,:));
        c2 = cross(q(1,3,:),q(1,4,:));
        alphak = sqrt(c1(1,1,1)^2 + c1(1,1,2)^2 + c1(1,1,3)^2)/...
            (sqrt(c1(1,1,1)^2+c1(1,1,2)^2+c1(1,1,3)^2) + sqrt(c2(1,1,1)^2 + c2(1,1,2)^2 + c2(1,1,3)^2) );
        alphak(isnan(alphak)) = 0;
        Vk = (1-alphak) .*q(1,2,:) + alphak .*q(1,3,:);
        Tu0l = Vk / sqrt((Vk(:,:,1))^2 + (Vk(:,:,2))^2 + (Vk(:,:,3))^2);
    end

%% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
    function Tukl = ComputeTukl (k, l)
        if ( l < 2 )
            q(1,3:4,:) = diff(Q(k+1,l+1:l+3,:),1,2);
            q(1,2,:) = 2 * q(1,3,:) - q(1,4,:);
            q(1,1,:) = 2 * q(1,2,:) - q(1,3,:);
        elseif ( l> m-2 )
            q(1,1:2,:) = diff(Q(k+1,l-2:l,:),1,2);
            q(1,3,:) = 2 * q(1,2,:) - q(1,1,:);
            q(1,4,:) = 2 * q(1,3,:) - q(1,2,:);
        else
            q = diff(Q(k+1,l-1:l+3,:),1,2);
        end
        c1 = cross(q(1,1,:),q(1,2,:));
        c2 = cross(q(1,3,:),q(1,4,:));
        alphak = sqrt(c1(1,1,1)^2 + c1(1,1,2)^2 + c1(1,1,3)^2)/...
            (sqrt(c1(1,1,1)^2+c1(1,1,2)^2+c1(1,1,3)^2) + sqrt(c2(1,1,1)^2 + c2(1,1,2)^2 + c2(1,1,3)^2) );
        alphak(isnan(alphak)) = 0;
        Vk = (1-alphak) .*q(1,2,:) + alphak .*q(1,3,:);
        Tukl = Vk / sqrt((Vk(:,:,1))^2 + (Vk(:,:,2))^2 + (Vk(:,:,3))^2);
    end

%% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)
    function Tvk0 = ComputeTvk0 (k)
        if ( k < 2 )
            q(3:4,1,:) = diff(Q(k+1:k+3,1,:),1,1);
            q(2,1,:) = 2 * q(3,1,:) - q(4,1,:);
            q(1,1,:) = 2 * q(2,1,:) - q(3,1,:);
        elseif ( k > n-2 )
            q(1:2,1,:) = diff(Q(k-2:k,1,:),1,1);
            q(3,1,:) = 2 * q(2,1,:) - q(1,1,:);
            q(4,1,:) = 2 * q(3,1,:) - q(2,1,:);
        else
            q = diff(Q(k-1:k+3,1,:),1,1);
        end
        c1 = cross(q(1,1,:),q(2,1,:));
        c2 = cross(q(3,1,:),q(4,1,:));
        alphak = sqrt(c1(1,1,1)^2 + c1(1,1,2)^2 + c1(1,1,3)^2)/...
            (sqrt(c1(1,1,1)^2+c1(1,1,2)^2+c1(1,1,3)^2) + sqrt(c2(1,1,1)^2 + c2(1,1,2)^2 + c2(1,1,3)^2) );
        alphak(isnan(alphak)) = 0;
        Vv = (1-alphak) .*q(2,1,:) + alphak .*q(3,1,:);
        Tvk0 = Vv / sqrt((Vv(:,:,1))^2 + (Vv(:,:,2))^2 + (Vv(:,:,3))^2);
    end

%% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
    function Tvkl = ComputeTvkl (k, l)
        if ( k < 2 )
            q(3:4,1,:) = diff(Q(k+1:k+3,l+1,:),1,1);
            q(2,1,:) = 2 * q(3,1,:) - q(4,1,:);
            q(1,1,:) = 2 * q(2,1,:) - q(3,1,:);
        elseif ( k > n-2 )
            q(1:2,1,:) = diff(Q(k-2:k,l+1,:),1,1);
            q(3,1,:) = 2 * q(2,1,:) - q(1,1,:);
            q(4,1,:) = 2 * q(3,1,:) - q(2,1,:);
        else
            q = diff(Q(k-1:k+3,l+1,:),1,1);
        end
        c1 = cross(q(1,1,:),q(2,1,:));
        c2 = cross(q(3,1,:),q(4,1,:));
        alphak = sqrt(c1(1,1,1)^2 + c1(1,1,2)^2 + c1(1,1,3)^2)/...
            (sqrt(c1(1,1,1)^2+c1(1,1,2)^2+c1(1,1,3)^2) + sqrt(c2(1,1,1)^2 + c2(1,1,2)^2 + c2(1,1,3)^2) );
        alphak(isnan(alphak)) = 0;
        Vv = (1-alphak) .*q(2,1,:) + alphak .*q(3,1,:);
        Tvkl = Vv / sqrt((Vv(:,:,1))^2 + (Vv(:,:,2))^2 + (Vv(:,:,3))^2);
    end


%% TODO: Load the U knot vector
%%       Load the V knot vector
    function [U, V] = LoadUV ()
        U = 0;
        V = 0;
        %Load the U knot vector (9.56)
        e = 8+2*(n-1);
        w = 8+2*(m-1);
        U = zeros(1,e);
        V = zeros(1,w);
        
        for j=1:4
            U(1,j) = 0;
        end
        for j=e-4:e
            U(1,j) = 1;
        end
        for j=4:e-4
            a = floor((j-1)/2);
            U(1,j)=ub(a);
        end
        
        %Load the V knot vector (9.56)
        for j=1:4
            V(1,j) = 0;
        end
        for j=w-4:w
            V(1,j) = 1;
        end
        for j=4:w-4
            a = floor((j-1)/2);
            V(1,j) = vb(a);
        end
    end




%% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k, l, 2)
    function Duvkl = ComputeDuvkl (k, l)

        if (l == 0)
            D_u(1,1:3,:) = r(l+1) .* td(k+1,1:3,:,1) ;
            Dvl = diff(vb(1:3));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            Delta_D_u = diff(D_u(1,1:3,:),1,2);
            dl = Delta_D_u ./Dvl';
            dl(:,2) = (1-bl) * dl(:,1) + bl * dl(:,2);
            c3 = 2;
            c4 = -1;
        elseif (l == m)
            D_u(1,1:3,:) = r(l+1) .* td(k+1,m-1:m+1,:,1) ;
            Dvl = diff(vb(m-1:m+1));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            Delta_D_u = diff(D_u(1,1:3,:),1,2);
            dl = Delta_D_u ./Dvl';
            dl(:,2) = (1-bl) * dl(:,1) + bl * dl(:,2);
            c3 = -1;
            c4 = 2;
        else
            D_u(1,1:3,:) = r(l+1) .* td(k+1,l:l+2,:,1) ;
            Dvl = diff(vb(l:l+2));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            Delta_D_u = diff(D_u(1,1:3,:),1,2);
            dl = Delta_D_u./Dvl';
            c3 = ( 1-bl);
            c4 = bl;
        end
        if (k == 0)
            D_v(1:3,1,:) = s(k+1) .* td(1:3,l+1,:,2) ;
            Duk = diff (ub(1:3));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            Delta_Dv = diff(D_v(1:3,1,:),1,1) ;
            dk  = Delta_Dv ./ Duk;
            dk(2,:)  = (1 - ak) * dk(1,:) + ak * dk(2,:);
            c1  =  2;
            c2  = -1;
            
        elseif (k == n)
            D_v(1:3,1,:) = s(k+1) .* td(n-1:n+1,l+1,:,2) ;
            Duk = diff (ub(n-1:n+1));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            Delta_Dv = diff(D_v(1:3,1,:),1,1) ;
            dk  = Delta_Dv ./ Duk;
            dk(2,:)  = (1 - ak) * dk(1,:) + ak * dk(2,:);
            c1  = -1;
            c2  =  2;
            
            
            
        else
            D_v(1:3,1,:) = s(k+1) .* td(k:k+2,l+1,:,2) ;
            Duk = diff (ub(k:k+2));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            Delta_Dv = diff(D_v(1:3,1,:),1,1);
            dk  = Delta_Dv ./ Duk;
            c1  = (1 - ak);
            c2  = ak;
            
        end
        
        
        d_vu = c1 * dk(1,:,:) + c2 * dk(2,:,:);
        d_uv = c3 * dl(:,1,:) + c4 * dl(:,2,:) ;
        Duvkl = (c2 * d_uv + c4 * d_vu) /(c2+c4);

    end

%% TODO: Compute the four inner control points of the (k,l)th
%% Bezier patch and load them into P
    function P = ComputeInnerControlPoints (P)
        gamma = (delta_uk(2:end)' * delta_vl(2:end))/9;
        for k=0:n-1
            for l=0:m-1
                
                P(3*k+2,3*l+2,:) = gamma * td(k+1,l+1,:,3) + P(k+1,l+2,:) + P(k+2,l+1,:) - P(k+1,l+1,:);
                P(3*k+3,3*l+2,:) = -gamma * td(k+2,l+1,:,3) + P(k+4,l+2,:) - P(k+4,l+1,:) + P(k+3,l+1,:);
                P(3*k+2,3*l+3,:) = -gamma * td(k+1,l+2,:,3) + P(k+2,l+4,:) - P(k+1,l+4,:) + P(k+1,l+3,:);
                P(3*k+3,3*l+3,:) = +gamma * td(k+2,l+2,:,3) + P(k+3,l+4,:) - P(k+4,l+3,:) + P(k+4,l+4,:);
                
            end
        end
    end

%% TODO: Load the NURBS control points by discarding Bezier points
%% along inner rows and columns, Figure 9.32c
    function P = DiscardInnerRowsandColumns (P)
        %Discarding rows
        for k=n-2:-1:0
            P(3*k+4,:,:)=[];
        end
        %Discarding columns
        for l=m-2:-1:0
            P(:,3*l+4,:)=[];
        end
    end

end
