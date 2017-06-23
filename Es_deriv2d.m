clc
clear all
close all

out = A93_94_95_test_data (10, 10 , 17, 25);
%Put ii=1 per prima figura,ii=2 per seconda
ii=2;
Q(:,:,1) = out{ii,1,1};
Q(:,:,2) = out{ii,1,2};
Q(:,:,3) = out{ii,1,3};
[n,m,o] = size(Q);

[uk,vl] = Surfmeshpar2(n,m,Q);
[Dk,Dl] = deriv_2d_for (uk,vl,Q) ;
figure
quiver3 (Q(:, :, 1), Q(:, :, 2), Q(:, :, 3), Dk(:, :, 1), Dk(:, :, 2), Dk(:, :, 3))
view (2)
hold on 
quiver3 (Q(:, :, 1), Q(:, :, 2), Q(:, :, 3), Dl(:, :, 1), Dl(:, :, 2), Dl(:, :, 3))
