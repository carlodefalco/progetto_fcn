function ders=DersBasisFunction(i,u,p,n,U)

ndu=ones(p+1);
left=zeros(1,p+1);
right=zeros(1,p+1);
for j=1:p
    left(j+1)=u-U(i+2-j);
    right(j+1)=U(i+j+1)-u;
    s=0;
    for r=0:j-1
        ndu(j+1,r+1)=right(r+2)+left(j-r+1);
        temp=ndu(r+1,j)/ndu(j+1,r+1);
        ndu(r+1,j+1)=s+right(r+2)*temp;
        s=left(j-r+1)*temp;
    end
    ndu(j+1,j+1)=s
end

ders=ones(p+1);


for j= 0:p
    ders(1,j+1)=ndu(j+1,p+1);
    for r=0:p
        s1=0;
        s2=1;
        a=ones(n+1);
        for k=1:n
            d=0;
            rk=r-k;
            pk=p-k;
            if r>=k
                a(s2+1,1)=a(s1+1,1)/ndu(pk+2,rk+1);
                d=a(s2+1,1)*ndu(rk+1,pk+1);
            end
                if rk>=-1
                    j1=1;
                else 
                    j1=-rk;
                end
                    if r-1<=pk
                        j2=k-1;
                    else
                        j2=p-r;
                    end
                        for j=j1:j2
                            a(s2+1,j+1)=(a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1);
                            d=d+a(s2+1,j+1)*ndu(rk+j+1,pk+1);
                        end
                            if r<=pk
                                a(s2+1,k+1)=-a(s1+1,k)/ndu(pk+2,r+1);
                                d=d+a(s2+1,k+1)*ndu(r+1,pk+1);
                            end
                            ders(k+1,r+1)=d;
                            j=s1;
                            s1=s2; 
                            s2=j;
                            

        end
    end

end

r=p;
for k=1:n
    for j=0:p
        ders(k+1,j+1)=ders(k+1,j+1)*r;
          r=r*(p-k);
        
    end
end

                           
                
            
        


                        
