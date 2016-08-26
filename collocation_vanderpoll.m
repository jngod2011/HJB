% soling vanderpoll oscillator HJB with cut points
clear
close all
clc


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


  
Q=diag([20,20,20]);
R=diag([1,1,1]);
% 
% Amat=[0,1;-1,1];
% Bmat=[0;1];
% [K,S,e] = lqr(Amat,Bmat,Q,R)



S=3*eye(d);
K=[2,2,2;2,2,2;2,2,2];


%%


dom=2;

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


% X=X+randn(size(X))
W=diag(w);

N=length(w);


%%
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
%%
V0=zeros(length(w),1);
for i=1:1:length(w)
    V0(i)=X(i,:)*S*X(i,:)';
end



%%
for ep=1
    

% % the lqr initial cost at these  collocation points


z=cell(length(w),1);
F=cell(length(w),1);
D=cell(length(w),1);

for i=1:1:length(w)
    xi=X(i,:)';
    D{i}=-1/4*B'*Jac(xi)*g(xi)*inv(R)*g(xi)'*Jac(xi)'*B;
    F{i}=B'*Jac(xi)*f(xi);
    z{i}=xi'*Q*xi;
end


options = optimset('Display','iter','Jacobian','on','Algorithm','levenberg-marquardt','MaxFunEvals',1e7,'MaxIter',1e7);
V0 = fsolve(@(V)quad_equal_const(V,D,F,z),V0,options)

% [VV,V0]
 [x,resnorm]=lsqnonlin(@(V)quad_equal_const(V,D,F,z),V0,zeros(1,N),[],options)
 
 keyboard
end

