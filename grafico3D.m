function [CC]=grafico3D(n,p,U,P)

t=linspace(0,1,100);
xx=linspace(-4,3,100);
CC=zeros(3,numel(t));

for ii=1:numel(t)
    span=findspan(n,p,t(ii),U);
    N=basisfun(span,t(ii),p,U);
    C=zeros(1,numel(P(1,:)));
    for i=0:p
    C=C+N(1,i+1)*P(span-p+i+1,:);
    end
    CC(:,ii)=C;
end
plot3(CC(1,:),CC(2,:),CC(3,:))

end