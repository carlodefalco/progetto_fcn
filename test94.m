clear all
clc
close all
w=20;
out = A93_94_95_test_data (w, w, 17, 25);

%metti ii=1 per cono ...
ii=1;
Q(:,:,1)=out{ii,1,1};
Q(:,:,2)=out{ii,1,2};
Q(:,:,3)=out{ii,1,3};
[n m o]=size(Q);
n=n-1;
m=m-1;
p=3;
q=3;

[uk,vl]=Surfmeshpar(n,m,Q);

[U,V,P]=Globalsurfinterp(n,m,Q,p,q);

knots = {U V} ;
cntrl(1,:,:)=P(:,:,1);
cntrl(2,:,:)=P(:,:,2);
cntrl(3,:,:)=P(:,:,3);
nrb = nrbmak(cntrl,knots);
subd= [40 40 ];
nrbplot(nrb,subd);
nrbplot(nrb,[w w]);

