uk = linspace (0, 2 * pi, 21);
Qk = [uk; sin(uk); zeros(size (uk))];
Dk = deriv_1d (uk, Qk);
n = size(uk,2) ;
dk = zeros(3,n-1);
Dk = zeros(size(Qk)) ;


for k=1:n-1
    Duk(k) = uk(k+1)-uk(k);
    qk(:,k) = Qk(:,k+1)-Qk(:,k) ;
    dk(:,k) = qk(:,k)./Duk(k) ;
end

for k=1:n-2
    ak(k) = Duk(k)/(Duk(k)+Duk(k+1)) ;
    Dk(:,k+1) = (1-ak(k))*dk(:,k) + ak(k)*dk(:,k+1) ;
end


for k=1:size(Qk,1)
    Dk(k,1)   = 2 * dk(k, 1)   - Dk(k, 2);
    Dk(k, end) = 2 * dk(k, end) - Dk(k, end-1);
end

figure
quiver3 (Qk(1, :), Qk(2, :), Qk(3, :), Dk(1, :), Dk(2, :), Dk(3, :))
view (2)