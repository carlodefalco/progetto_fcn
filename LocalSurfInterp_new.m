function [U, V, P] = LocalSurfInterp_formaprof (n, m, Q)
%% get ub, r and u direction tangents
td=zeros(n+1,m+1,3);
total = 0;
ub = zeros (n+1, 1);
for l=0:m
    %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    Tu=ComputeTu0l(n,m,Q,l)
    td(1,l+1,1)=Tu(l+1);
    r(l+1)=0;
    for k=1:n
        %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
        Tu=ComputeTukl(n,m,Q,l,k)
        td(2:end,l+1,1)=Tu(k+1,l+1);
        d =sqrt((Q(l+1, k+1, 1) - Q(l+1, k, 1))^2+(Q(l+1, k+1, 2) - Q(l+1, k, 2))^2+(Q(l+1, k+1, 2) - Q(l+1, k+1, 2))^2);
        ub(k+1) = ub(k+1) + d;
        r(l+1)  = r(l+1)  + d;
    end
    total = total + r(l+1) ;
end

for k=1:n
    ub(k+1) = ub(k) + ub(k+1) / total ;
    ub(n+1) = 1;
end

%% get vb, s and v direction tangents
total = 0;
vb = zeros (m+1, 1);
for k = 0:n
    %% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)
    Tv = ComputeTvk0 (n,m,Q,k);
    td(k+1,1,2)=Tv(k+1);
    s(k+1) = 0;
    for l = 1:m
        %% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
        Tv = ComputeTvkl (n,m,Q,k,l);
        td(k+1,2:end,2)=Tv(k+1,l+1)
        d =sqrt((Q(l+1, k+1, 1) - Q(l, k+1, 1))^2+(Q(l+1, k+1, 2) - Q(l, k+1, 2))^2+(Q(l+1, k+1, 2) - Q(l, k+1, 2))^2);
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
 
  [U, V] = LoadUV (ub,vb,n,m);
  
 %%       Compute all Bezier control points along each
 %%       row and column of data points.
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
        
        a=r(k+1)*delta_uk(i+1)/3;
        P(3*i+2,3*k+1,:)=Q(i+1,k+1,:)+a*td(i+1,k+1,1);
        P(3*i+3,3*k+1,:)=Q(i+2,k+1,:)-a*td(i+2,k+1,1);
    end
end
%Ciclo for per le righe di Bezier control points

for i=0:m
    for k=0:n-1
        b=s(k+1)*delta_vl(k+1)/3;
        P(3*i+1,3*k+2,:)=Q(i+1,k+1,:)+b*td(i+1,k+1,2);
        P(3*i+1,3*k+3,:)=Q(i+1,k+2,:)-b*td(i+1,k+2,2);
    end
end



  for k = 0:n
    for l = 0:m
      %% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k+1, l+1, 2+1)
      Duvkl = ComputeDuvkl(n,m,Q,r,s,td,delta_uk,delta_vl);
      td(k+1, l+1, 3)=Duvkl(k+1,l+1);
    end
  end

  
      %% TODO: Compute the four inner control points of the (k,l)th
      %% Bezier patch and load them into P

gamma=(delta_uk*delta_vl.')/9;

for k=0:n-1
    for i=0:m-1
        
        
        P(3*i+2,3*k+2,:)=gamma*Duvkl(i+1,k+1)+P(i+1,k+2,:)+P(i+2,k+1,:)-P(i+1,k+1,:);
        P(3*i+3,3*k+2,:)=-gamma*Duvkl(i+2,k+1)+P(i+4,k+2,:)-P(i+4,k+1,:)+P(i+3,k+1,:);
        P(3*i+2,3*k+3,:)=-gamma*Duvkl(i+1,k+2)+P(i+2,k+4,:)-P(i+1,k+4,:)+P(i+1,k+3,:);
        P(3*i+3,3*k+3,:)=+gamma*Duvkl(i+2,k+2)+P(i+3,k+4,:)-P(i+4,k+3,:)+P(i+4,k+4,:);
        
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

end
