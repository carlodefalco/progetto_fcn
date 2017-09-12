clc
clear all
close all

out = A93_94_95_test_data (7 , 5 , 17, 25);
%Put ii=1 per prima figura,ii=2 per seconda
ii = 2 ;
Q(:,:,1)=out{ii,1,1};
Q(:,:,2)=out{ii,1,2};
Q(:,:,3)=out{ii,1,3};
[n,m,o]=size(Q);
n=n-1;
m=m-1;

[U,V,P,td]=LocalSurfInterp(n,m,Q);

knots = {U V} ;
cntrl(1,:,:)=P(:,:,1);
cntrl(2,:,:)=P(:,:,2);
cntrl(3,:,:)=P(:,:,3);
nrb = nrbmak(cntrl,knots);
figure
nrbplot(nrb,[n m]);
% 
% hold on
% plot3(P(:,:,1),P(:,:,2),P(:,:,3),'x')
% plot3(Q(:,:,1),Q(:,:,2),Q(:,:,3),'o')
% hold off
% 
% hold on
% mesh(P(:,:,1),P(:,:,2),P(:,:,3))
% % mesh(Q(:,:,1),Q(:,:,2),Q(:,:,3))
% hold off
