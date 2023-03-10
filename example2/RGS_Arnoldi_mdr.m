function [Q,R,Z,S]=RGS_Arnoldi_mdr(A,b,theta,m,M1,M2,Vn,Sn)
n=length(b);
Q=zeros(n,m+1);
R=zeros(m+1,m);
Z=zeros(n,m);

w=b;
p=theta*w;
newq=w;
news=p;
h=norm(news);
S(:,1)=news/h;
Q(:,1)=newq/h;

for i=1:m
    if isempty(M1)==0
        Z(:,i)=M1\Q(:,i);
        if isempty(M2)==0
            Z(:,i)=M2\Z(:,i);
        end
    else
       Z(:,i)=Q(:,i); 
    end
    w=A*Z(:,i);
    w=w-Vn*(Sn'*(theta*w));
    p=theta*w;
    R(1:i,i)=S(:,1:i)\p;
    newq=w-Q(:,1:i)*R(1:i,i);
    news=theta*newq;
    R(i+1,i)=norm(news,2);
    S(:,i+1)=news/R(i+1,i);
    Q(:,i+1)=newq/R(i+1,i);
end
end