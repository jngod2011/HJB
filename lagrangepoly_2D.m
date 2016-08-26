% 2D lagrange polynomials


n=2;
N=3;

x1=-2;
x2=2;
x3=0;
y1=-1;
y2=2;
y3=0;

S=[4,1;1,4];
fx1y1=[x1,y1]*S*[x1;y1];
fx2y1=[x2,y1]*S*[x2;y1];
fx1y2=[x1,y2]*S*[x1;y2];
fx2y2=[x2,y2]*S*[x2;y2];
fx3y3=[x3,y3]*S*[x3;y3];
% [x1,y1]
phi1=@(x,y,x1,x2,x3,y1,y2,y3)    (x-x2)*(x-x3)/((x1-x2)*(x1-x3)) * (y-y2)*(y-y3)/((y1-y2)*(y1-y3))    ;
% [x2,y1]
phi2=@(x,y,x1,x2,x3,y1,y2,y3)    (x-x1)*(x-x3)/((x2-x1)*(x2-x3)) * (y-y2)*(y-y3)/((y1-y2)*(y1-y3))    ;
% [x1,y2]
phi3=@(x,y,x1,x2,x3,y1,y2,y3)    (x-x2)*(x-x3)/((x1-x2)*(x1-x3)) * (y-y1)*(y-y3)/((y2-y1)*(y2-y3))    ;
% [x2,y2]
phi4=@(x,y,x1,x2,x3,y1,y2,y3)    (x-x1)*(x-x3)/((x2-x1)*(x2-x3)) * (y-y1)*(y-y3)/((y2-y1)*(y2-y3))    ;

% [x3,y3]
phi5=@(x,y,x1,x2,x3,y1,y2,y3)    (x-x1)*(x-x2)/((x3-x1)*(x3-x2)) * (y-y1)*(y-y2)/((y3-y1)*(y3-y2))    ;

[xx,yy]=meshgrid(-4:0.1:4);
p=zeros(size(xx));
pquad=p;
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        p(i,j)=fx1y1*phi1(xx(i,j),yy(i,j),x1,x2,x3,y1,y2,y3)+ fx2y1*phi2(xx(i,j),yy(i,j),x1,x2,x3,y1,y2,y3) +fx1y2*phi3(xx(i,j),yy(i,j),x1,x2,x3,y1,y2,y3)+fx2y2*phi4(xx(i,j),yy(i,j),x1,x2,x3,y1,y2,y3)...
            +fx3y3*phi5(xx(i,j),yy(i,j),x1,x2,x3,y1,y2,y3);
        pquad(i,)
    end
end


mesh(xx,yy,p)
hold on
X=[x1,y1,fx1y1; x2,y1,fx2y1; x1,y2,fx1y2; x2,y2,fx2y2];
plot3(X(:,1),X(:,2),X(:,3),'ro','MarkerSize',8,'linewidth',2)




% z= ax+ by+c 
A=[xy1(1),xy1(2),1;
   xy2(1),xy2(2),1;
   xy3(1),xy3(2),1];
syms x y
M1=A;
M1(1,1)=1;
M1(1,2)=1;
M1=M1.*[x,y,1;1,1,1;1,1,1];
M1=det(M1);

M2=A;
M2(2,1)=1;
M2(2,2)=1;
M2=M2.*[1,1,1;x,y,1;1,1,1];
M2=det(M2);

M3=A;
M3(3,1)=1;
M3(3,2)=1;
M3=M3.*[1,1,1;1,1,1;x,y,1];
M3=det(M3);

M=det(A);

f=f1*M1/M+f2*M2/M+f3*M3/M;


[xx,yy]=meshgrid(0:0.1:4);
fp=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        fp(i,j)=subs( subs(f,'x',xx(i,j)) ,'y',yy(i,j));
    end
end

mesh(xx,yy,fp)
hold on
X=[xy1,f1;xy2,f2;xy3,f3];
plot3(X(:,1),X(:,2),X(:,3),'ro','MarkerSize',8,'linewidth',2)


