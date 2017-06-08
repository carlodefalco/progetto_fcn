out = A92_test_data (n, 100)


%Qui metti 1 per prima figura, 2 per seconda etc etc...
ii=4;
Q=[out{ii,1}(1,:);out{ii,1}(2,:);out{ii,1}(3,:)]';



x=Q(:,1);
y=Q(:,2);
z=Q(:,3);
n=numel(x)-1;
r=1;
p=3;
[U,m,P]=alg1_andri(n,p,Q,r);

for i=1:n
d(i)=abs([x(i+1)-x(i)]^2+[y(i+1)-y(i)]^2+[z(i+1)-z(i)]^2)^(1/2);
end
d_tot=sum(d);

%Calcolo dei parametri u segnato
u=zeros(n+1,1)';
u(1)=0;
u(numel(x))=1;

for i=1:(n-1)
u(i+1)=u(i)+d(i)/d_tot;
end
U=zeros(1,n+7)
p=3;
for i=1:p+1
end
for i=n+4:n+7
U(1,i)=1
end
for i=2:n
U(1,i+3)=u(i)
end

%u1 e d sono due input per le funzioni che seguono
u1=1/4;
d=4;
du=min(d,p)+1;

%CK=zeros(d,2);
for k=p+1:d
    CK(k+1,3)=0;
    span=findspan(n,p,u1,U);
    nders=DersBasisFunction(span,u1,p,du,U);
    for k=0:du
        CK(k+1,3)=0;
        for j=0:p
          CK(k+1,:)=CK(k+1,:)+nders(k+1,j+1)*P(span-p+j+1,:);
        end
    end
end

%Calcolo dei primi ed ultimi due punti di P
P=zeros(n+3,3);
P(1,:)=Q(1,:);
P(n+3,:)=Q(n+1,:);
P(2,:)=P(1,:)+U(5)/3*CK(1,:);
P(n+2,:)=P(n+3,:)-(1-U(n+3))/3*CK(d+1,:);

% %Ora inizia il vero algoritmo

R=ones(n+1,3);
for i=3:n-1
    R(i+1,:)=Q(i,:);
    abc=basisfuncs(4,U(5),3,U);
    den=abc(2);
    P(3,:)=(Q(2,:)-abc(1)*P(2,:))/den;
    for i=3:n-1 
        dd(i+1)=abc(3)/den;
        abc=basisfuncs(i+2,U(i+3),3,U);
        den=abc(2)-abc(1)*dd(i+1);
        P(i+1,:)=(R(i+1,:)-abc(1)*P(i,:))/den;
    end
    dd(n+1)=abc(3)/den;
    abc=basisfuncs(n+2,U(n+3),3,U);
    den=abc(2)-abc(1)*dd(n+1);
    P(n+1,:)=(Q(n,:)-abc(3)*P(n+2,:)-abc(1)*P(n,:))/den;
    for i=n-1:-1:2
        P(i+1,:)=P(i+1,:)-dd(i+2)*P(i+2,:);
    end
end


%Plot della funzione tramite algoritmo da noi scritto

figure
plot3(Q(:,1),Q(:,2),Q(:,3),'o')
hold on, box on


[CC]=grafico3D(n,p,U,P);
    
