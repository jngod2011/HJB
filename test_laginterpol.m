
clc
clear

S=[4,0;0,4];

dom=2;

Ncol_gh=3;

d=2;

bdd_low=-dom*ones(1,d);
bdd_up=dom*ones(1,d);

[X,w] = GLeg_pts(Ncol_gh*ones(1,d),bdd_low,bdd_up);


Pbasis=InterpPoly(X);
ss=GenMfile_MatrixOfPolys(Pbasis,'','');
phi = str2func(strcat('@(x)',ss));

Vbars=zeros(length(w),1);
for i=1:1:length(w)
    Vbars(i)=X(i,:)*S*X(i,:)';
end



[xx,yy]=meshgrid(linspace(-dom,dom,20));
Vlag=zeros(size(xx));
Vtr=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        Vlag(i,j)=Vbars'*phi([xx(i,j),yy(i,j)]);
        Vtr(i,j)=[xx(i,j),yy(i,j)]*S*[xx(i,j);yy(i,j)];
    end
end

mesh(xx,yy,Vtr)
alpha 0.4
hold on

mesh(xx,yy,Vlag)
% alpha 0.4
plot3(X(:,1),X(:,2),Vbars,'ro','linewidth',2,'MarkerSize',7)


