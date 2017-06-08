function [uk,vl]=Surfmeshpar1(n,m,Q)

x=[Q(:,:,1)];
y=[Q(:,:,2)];
z=[Q(:,:,3)];

num=m+1; 
uk=zeros(1,n+1);
uk(1,1)=0;
uk(1,n+1)=1;


for k=1:n-1
    uk(1,k+1)=0;
end
for l=0:m
        total=0; %Total chord length of the row
        for k=1:n
            %cds=pdist(Q(k+1,l+1),Q(k,l+1),'euclidean')
            cds(k)=sqrt((x(k+1,l+1)-x(k,l+1))^2+(y(k+1,l+1)-y(k,l+1))^2+(z(k+1,l+1)-z(k,l+1))^2);     %cds(k+1)=distance3D(Q(k+1,l),Q(k,l))
            total=total+cds(k);
        end
        if total==0
            num=num-1;
        else
            d=0;
            for k=1:n-1
                d=d+cds(k);
                uk(1,k+1)=uk(1,k+1)+d/total;
            end
        end
    end
    if num==0 
        display('Error (num=0)')
    end
    for k=1:n-1
        uk(1,k+1)=uk(1,k+1)/num;
    end
%Now compute the vl

num=n+1;
vl=zeros(1,m+1);
vl(1,1)=0;
vl(1,m+1)=1;

for l=1:m-1
    vl(1,l+1)=0;
end
    for k=0:n
        total=0;
        for l=1:m
            %cds=pdist(Q(k+1,l+1),Q(k,l+1),'euclidean')
            cds(l)=sqrt((x(k+1,l+1)-x(k+1,l))^2+(y(k+1,l+1)-y(k+1,l))^2+(z(k+1,l+1)-z(k+1,l))^2);     %cds(k+1)=distance3D(Q(k+1,l),Q(k,l))
            total=total+cds(l);
        end
            if total==0
            num=num-1;
        else
            d=0;
            for l=1:m-1
                d=d+cds(l+1);
                vl(1,l+1)=vl(1,l+1)+d/total;
            end
        end
    end
    if num==0 
        display('Error (num=0)')
    end
    for l=1:m-1
        vl(1,l+1)=vl(1,l+1)/num;
    end

end


    

