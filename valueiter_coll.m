% value iteration by collocation
clc
clear all
close all
%% Given system is
Q=diag([20,20,20]);
R=diag([1,1,1]);


d=3;

% % ep=1;
% g=@(x)[0;1];
% f=@(x)[x(2);-x(1)+x(2)-ep*x(1)^2*x(2)];
%f=@(x)[x(2);x(1)-x(1)^3];

I1=10;
I2=15;
I3=20;
g=@(x)eye(d);
f=@(x)[(I2-I3)/I1 * x(2)*x(3);...
      (I3-I1)/I2 * x(3)*x(1);...
      (I1-I2)/I3 * x(1)*x(2)];


dom=2;

K=[2,0,0;0,2,0;0,0,2];
u=@(x)-f(x)-K*x(:);
[xx,yy]=meshgrid(linspace(-dom,dom,20));

%%
% [xx,yy]=meshgrid(linspace(-dom,dom,20));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%     [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,100],[xx(i,j);yy(i,j);1]);
%     plot(x(:,1),x(:,2))
%     hold on
%     end
% end


%% LQR solution


% Amat=[0,1;-1,1];
% Bmat=[0;1];
% [K,S,e] = lqr(Amat,Bmat,Q,R);

S=3*eye(d);
K=[2,2,2;2,2,2;2,2,2];
%% the points that will be used for collocation


Ncol_gh=3;

bdd_low=-dom*ones(1,d);
bdd_up=dom*ones(1,d);
X=[];
w=[];
X2=[];
w2=[];
% [X,w]=uniform_sigma_pts(bdd_low,bdd_up,6);
% [X2,w2]=uniform_sigma_pts(bdd_low/2,bdd_up/2,4);
[X,w] = GLeg_pts(Ncol_gh*ones(1,d),bdd_low,bdd_up);
[X2,w2] = GLeg_pts((Ncol_gh-1)*ones(1,d),bdd_low/2,bdd_up/2);

X=[X;X2];
w=[w;w2];

X=X+0.01*randn(size(X));

% X=X+randn(size(X))
W=diag(w);

N=length(w);

% % the lqr initial cost at these  collocation points

V0=zeros(length(w),1);
for i=1:1:length(w)
    V0(i)=X(i,:)*S*X(i,:)';
end


%% generating all the polynomia basis and the correponding matrices here

order_poly_approx=5;

Pbasis=Basis_polyND(d,order_poly_approx);
ss=GenMfile_MatrixOfPolys(Pbasis,'','');
phi = str2func(strcat('@(x)',ss));


Jac=Jacobian_VectorOfPolys(Pbasis);

ss=GenMfile_MatrixOfPolys(Jac,'','');
Jac = str2func(strcat('@(x)',ss));

m=length(Pbasis);

PPT=AddSubMultiply_MatrixOfPolys(Pbasis,Pbasis','multiply');
ss=GenMfile_MatrixOfPolys(PPT,'','');
PPT = str2func(strcat('@(x)',ss));

A=zeros(N,m);
for i=1:1:N
    A(i,:)=phi(X(i,:));
end




BmapType='interp';
% BmapType='2norm';

if strcmpi(BmapType,'2norm')
    % check for simple quadratic
    [XX,ww] = GLeg_pts(11*ones(1,d),bdd_low,bdd_up);
    M=0;
    for i=1:1:length(ww)
        M=M+ww(i)*PPT(XX(i,:));
    end
    B=inv(M)*A'*W;
    
elseif strcmpi(BmapType,'interp')
    
    if N<m
    B=A'*inv(A*A');
    else
    B=inv(A'*A)*A';     
    end
end

% keyboard
% u=@(x)-K*x(:);

c=randn(m,1).^2;



type='Ccoeff';
% type='Vcoeff';
% type='Galerkin';


for iter=1:1:50
    
    
    
    
    if strcmpi(type,'Vcoeff')
        D=[];
        F=[];
        for i=1:1:N
            xi=X(i,:)';
            Jxi=Jac(X(i,:));
            D=vertcat(D,     (f(xi)+g(xi)*u(xi))'*Jxi'*B   ) ;
            F=vertcat(F,     -xi'*Q*xi-u(xi)'*R*u(xi)      ) ;
        end
        ind=find(sum(abs(D),2)==0);
        D(ind,:)=[];
        F(ind,:)=[];
        
%                 keyboard
        %Vc=D\(F);
        %  c=Apinv*Vc;
        
        % derivarive constraint
        D=vertcat(D,Jac(zeros(1,d))'*B);
        F=vertcat(F,zeros(d,1));
        
        %     % 0 at the origin constraint
        D=vertcat(D,phi(zeros(1,d))'*B);
        F=vertcat(F,0);
        
        options = optimset('Display','iter','Jacobian','off','Algorithm','interior-point','MaxFunEvals',1e7,'MaxIter',1e7);
        Vb=fmincon(@(Vb)sum(Vb.^2),V0,[],[],D,F,zeros(1,N),[],[],options)   %
        c=B*Vb;  
        
    elseif strcmpi(type,'Ccoeff')
        
        D=[];
        F=[];
        Aineq=-A;
        for i=1:1:N
            xi=X(i,:)';
            Jxi=Jac(X(i,:));
%             Aineq=vertcat(Aineq,f(xi)'*Jxi');
            D=vertcat(D,     (f(xi)+g(xi)*u(xi))'*Jxi'   ) ;
            F=vertcat(F,     -xi'*Q*xi-u(xi)'*R*u(xi)      ) ;
        end
%         keyboard
        %Vc=D\(F);
        %  c=Apinv*Vc;
        
        % derivarive constraint
%                 D=vertcat(D,Jac([0,0])');
%                 F=vertcat(F,zeros(d,1));
        
        %     % 0 at the origin constraint
        D=vertcat(D,phi(zeros(1,d))');
        F=vertcat(F,0);
        
        %adding the inequality constraints
        
        
        
        
        options = optimset('Display','iter','Jacobian','off','Algorithm','interior-point','MaxFunEvals',1e7,'MaxIter',1e7);
        c=fmincon(@(c)sum(c.^2),c,Aineq,zeros(size(Aineq,1),1),D,F,[],[],[],options)
        
        
    elseif strcmpi(type,'Galerkin')
        
        D=[];
        F=[];
        for i=j:1:m
            for i=1:1:N
                xi=X(i,:)';
                Jxi=Jac(X(i,:));
                phi(xi')
                
                D=vertcat(D,     (f(xi)+g(xi)*u(xi))'*Jxi'   ) ;
                F=vertcat(F,     -xi'*Q*xi-u(xi)'*R*u(xi)      ) ;
            end
        end
        
        
    end
    
    VVlqr=zeros(size(xx));
    VVapprox=zeros(size(xx));
    err=zeros(size(xx));
    for i=1:1:size(xx,1)
        for j=1:1:size(xx,2)
            xi=[xx(i,j);yy(i,j);1];
            err(i,j)=HJBerror(Q,R,f,g,@(x)Jac(x),xi,c);
            VVlqr(i,j)=xi'*S*xi;
            VVapprox(i,j)=c'*phi(xi');
        end
    end
    figure(1)
    subplot(1,2,1)
    hold off
    %     mesh(xx,yy,VVlqr)
    %     alpha 0.4
    
    mesh(xx,yy,VVapprox)
    hold on
    alpha 0.4
    plot(X(:,1),X(:,2),'ro','linewidth',2,'MarkerSize',6)
    
    subplot(1,2,2)
    hold off
    mesh(xx,yy,err)
    hold on
    alpha 0.4
    plot(X(:,1),X(:,2),'ro','linewidth',2,'MarkerSize',6)
    
    pause(1)
    
    u=@(x)-0.5*inv(R)*g(x)'*Jac(x)'*c;
    
    
    
    
    
end


err=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        xi=[xx(i,j);yy(i,j);1];
        err(i,j)=HJBerror(Q,R,f,g,@(x)Jac(x),xi,c);
    end
end
figure(2)
mesh(xx,yy,err)


[xx,yy]=meshgrid(linspace(-dom,dom,10));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);1]);

        plot(x(:,1),x(:,2))
        hold on
    end
end
