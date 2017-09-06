clear all
clc
close all

out = A93_94_95_test_data (100, 100, 17, 25);
ii=1;
Q(:,:,1)=out{1,1,1};
Q(:,:,2)=out{1,1,2};
Q(:,:,3)=out{1,1,3};
[r s o]=size(Q);
r=r-1;
s=s-1;
n=15; %Number of control points
m=8;
p=3;
q=3;

x=Q(:,:,1);
y=Q(:,:,2);
z=Q(:,:,3);

 [U,V,P]=GlobSurfApprox(r,s,Q,p,q,n,m)
 figure
knots = {U V} ;
% cntrl = zeros [3 , n , m];
cntrl(1,:,:)=P(:,:,1);
cntrl(2,:,:)=P(:,:,2);
cntrl(3,:,:)=P(:,:,3);
nrb = nrbmak(cntrl,knots);
nrbplot(nrb,[20 20]);
