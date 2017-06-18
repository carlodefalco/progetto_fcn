function [Dk,Dl] = deriv_2d_for (uk, vl, Qk)

n = numel (uk) - 1 ;
m = numel(vl) -1 ;

for k = 0 : n
    for l = 0 : m
        if (l == 0)
            Dvl = diff(vl(1:3));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            ql = diff (Qk(1:3,:),1,1);
            dl = ql./Dvl';
            dl(2,:) = (1-bl) * dl(1,:) + bl * dl(2,:);
            c3 = 2;
            c4 = -1;
        elseif (l==m)
            Dvl = diff(vl(m-1:m+1));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            ql = diff (Qk(m-1:m+1,:),1,1);
            dl = ql./Dvl';
            dl(2,:) = (1-bl) * dl(1,:) + bl * dl(2,:);
            c3 = -1;
            c4 = 2;
        else
            Dvl = diff(vl(l:l+2));
            bl = Dvl(1) / (Dvl(1) + Dvl(2));
            ql = diff (Qk(l:l+2,:),1,1);
            dl = ql./Dvl';
            c3 = ( 1-bl);
            c4 = bl;
        end
        if (k == 0)
            
            Duk = diff (uk(1:3));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            qk  = diff (Qk(:, 1:3), 1, 2);
            dk  = qk ./ Duk;
            dk(:,2)  = (1 - ak) * dk(:,1) + ak * dk(:,2);
            c1  =  2;
            c2  = -1;
            
            
            
        elseif (k == n)
            
            Duk = diff (uk(n-1:n+1));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            qk  = diff (Qk(:, n-1:n+1), 1, 2);
            dk  = qk ./ Duk;
            dk(:,1) = (1 - ak) * dk(:,1) + ak * dk(:,2);
            c1  = -1;
            c2  =  2;
            
            
            
        else
            
            Duk = diff (uk(k:k+2));
            ak  = Duk(1) / (Duk(1) + Duk(2));
            qk  = diff (Qk(:, k:k+2), 1, 2);
            dk  = qk ./ Duk;
            c1  = (1 - ak);
            c2  = ak;
            
        end
        
        
        Dk(:, k+1) = c1 * dk(:,1) + c2 * dk(:,2);
        Dl(l+1,:) = c3 * dl(1,:) + c4 * dl(2,:) ;
        
    end
end
end
