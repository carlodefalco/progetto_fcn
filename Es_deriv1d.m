clc
clear all

uk = linspace (0, 2 * pi, 21);
Qk = [uk; sin(uk); zeros(size (uk))];
n = size(uk,2) ;
Dk = zeros(size(Qk)) ;

for k=0:n-3
        Duk(1) = uk(2)-uk(1);
        Duk(k+2) = uk(k+3)-uk(k+2);
        ak = Duk(k+1)/(Duk(k+1)+Duk(k+2));
        qk(:,1) = Qk(:,2)-Qk(:,1) ;
        qk(:,k+2) = Qk(:,k+3)-Qk(:,k+2) ;
        dk(:,1) = qk(:,1)./Duk(1);
        dk(:,k+2) = qk(:,k+2)./Duk(k+2);
        Dk(:,k+2) = (1-ak)*dk(:,k+1) + ak*dk(:,k+2) ;
end

Dk(:,1)   = 2 * dk(:, 1)   - Dk(:, 2);
Dk(:, end) = 2 * dk(:, end) - Dk(:, end-1);

figure
quiver3 (Qk(1, :), Qk(2, :), Qk(3, :), Dk(1, :), Dk(2, :), Dk(3, :))
view (2)
