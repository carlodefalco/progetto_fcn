clc
clear all
close all
out = A93_94_95_test_data (10 , 10, 17, 25);
%Put ii=1 per prima figura,ii=2 per seconda
ii=2;
Qx=out{ii,1,1};
Qy=out{ii,1,2};
Qz=out{ii,1,3};
[n,m]=size(Qx);
n=n-1;
m=m-1;
 
[Ux,Vx,Px]=LocalSurfInterp_formaprof(n,m,Qx);
[Uy,Vy,Py]=LocalSurfInterp_formaprof(n,m,Qy);
[Uz,Vz,Pz]=LocalSurfInterp_formaprof(n,m,Qz);
P(:,:,1)=Px;
P(:,:,2)=Py;
P(:,:,3)=Pz;

