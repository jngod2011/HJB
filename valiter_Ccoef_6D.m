% value iteration by collocation
clc
clear all
close all
%% Given system is
Q=diag([40,40,40,40,40,40]);
R=diag([1,1,1]);

figname='6DAttStab_lowCtrlCost_CUT_';
saveonoff=1;



d=6;

Ix=14;
Iy=10;
Iz=8;

% rodrigues parameters by crassidis
% g=@(x)[zeros(3);diag([1/Ix,1/Iy,1/Iz])];
% f=@(x)[0.5*(0.5*(1-[x(1),x(2),x(3)]*[x(1);x(2);x(3)])*eye(3)-[0,x(6),-x(5);-x(6),0,x(4);x(5),-x(4),0]+[x(1);x(2);x(3)]*[x(1),x(2),x(3)])*[x(4);x(5);x(6)];
%     x(5)*x(6)*(Iy-Iz)/Ix;
%     x(4)*x(6)*(Iz-Ix)/Iy;
%     x(4)*x(5)*(Ix-Iy)/Iz];
%
%
% dom=[-2,2;
%     -2,2;
%     -2,2;
%     -1,1;
%     -1,1;
%     -1,1];


% rodrigues parameters 2 by beard
g=@(x)[zeros(3);diag([1/Ix,1/Iy,1/Iz])];
f=@(x)[0.5*(eye(3)-[0,x(6),-x(5);-x(6),0,x(4);x(5),-x(4),0]+[x(1);x(2);x(3)]*[x(1),x(2),x(3)])*[x(4);x(5);x(6)];
    x(5)*x(6)*(Iy-Iz)/Ix;
    x(4)*x(6)*(Iz-Ix)/Iy;
    x(4)*x(5)*(Ix-Iy)/Iz];


dom=[-2,2;
    -2,2;
    -2,2;
    -1,1;
    -1,1;
    -1,1];



% euler angles
% g=@(x)[zeros(3);diag([1/Ix,1/Iy,1/Iz])];
% f=@(x)[x(4)+x(6)*cos(x(1))*tan(x(2))+x(5)*sin(x(1))*tan(x(2));
%     x(5)*cos(x(1))-x(6)*sin(x(1));
%     x(6)*cos(x(1))*sec(x(2))+x(5)*sec(x(2))*sin(x(1));
%     x(5)*x(6)*(Iy-Iz)/Ix;
%     x(4)*x(6)*(Iz-Ix)/Iy;
%     x(4)*x(5)*(Ix-Iy)/Iz];
%
%
% dom=[-60*pi/180,60*pi/180;
%     -60*pi/180,60*pi/180;
%     -60*pi/180,60*pi/180;
%     -1,1;
%     -1,1;
%     -1,1];


% get the optimized traj for this initial state
x0=[1.4735,0.6115,2.5521,0,0,0];


%%


K=[-15,0,0,-5,0,0;
    0,-15,0,0,-5,0;
    0,0,-15,0,0,-5];
u=@(x)K*x(:);
u0=@(x)K*x(:);

% close all
% figure
% [xx,yy]=meshgrid(linspace(dom(1,1),dom(1,2),10)*0.8,linspace(dom(2,1),dom(2,2),10)*0.8);
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);dom(3,2)*0.8;dom(4,2)*0.8;dom(5,2)*0.8;dom(6,2)*0.8]);
%         [i,j]
%         plot(x(:,1),x(:,2),'linewidth',2)
%         norm(x(end,:))
%         axis([-2,2,-2,2])
%         hold on
%         pause(0.2)
%     end
% end
% xlabel('p_1')
% ylabel('p_2')
% plot_prop_paper
% if saveonoff==1
%     saveas(gca,strcat(figname,'InitialTrajs'),'eps')
%     saveas(gca,strcat(figname,'InitialTrajs'),'png')
%     saveas(gca,strcat(figname,'InitialTrajs'),'fig')
% end
%% the points that will be used for collocation


Ncol_gh=3;

% bdd_low=-dom*ones(1,d);
% bdd_up=dom*ones(1,d);
X=[];
w=[];
X2=[];
w2=[];
X3=[];
w3=[];
X4=[];
w4=[];

% [X,w]=smolyak_sparse_grid_modf(zeros(d,1),diag([1.3,1.3,1.3,0.8,0.8,0.8]),d,5,'GH');
% mu=(dom(:,1)'+dom(:,2)')/2;
% h=-dom(:,1)'+dom(:,2)';
% for i=1:1:d
%     X(:,i)=(h(i)/2)*X(:,i)+mu(i);
% end

[X,w]=conjugate_dir_gausspts_till_8moment(zeros(d,1),diag([1.5,1.5,1.5,1,1,1]));
%  4.0249    4.0249    4.0249    3.2863    3.2863    3.2863
% [X,w] = GH_points(zeros(d,1),diag([5.5,5.5,5.5,3.5,3.5,3.5]),3);
% [X2,w2]=conjugate_dir_gausspts_till_6moment_scheme2(zeros(d,1),diag([0.6,0.6,0.6,0.2,0.2,0.2]));
max(abs(X),[],1)
% [X,w]=uniform_sigma_pts(dom(:,1)',dom(:,2)',6);

[X2,w2]=uniform_sigma_pts(dom(:,1)'/3,dom(:,2)'/3,6);
% [X3,w3] = GLeg_pts(Ncol_gh*ones(1,d),dom(:,1)'/2,dom(:,2)'/2);
% [X3,w3] = GLeg_pts((Ncol_gh)*ones(1,d),dom(:,1)'*0.5,dom(:,2)'*0.5);
% [X4,w4] = GLeg_pts(Ncol_gh*ones(1,d),dom(:,1)'*0.25,dom(:,2)'*0.25);

% X3=[general_conj_axis(d,d);general_conj_axis(d,2)]
% 
% for i=1:1:d
%     X3(X3(:,i)==-1,i)=dom(i,1);
%     X3(X3(:,i)==1,i)=dom(i,2);
% end
% w3=ones(size(X3,1),1)/size(X3,1);

X=[X;X2;X3;X4];
w=[w;w2;w3;w4];

W=diag(w);

N=length(w);


% Nrnd=200;
% Xrnd=(rand(Nrnd,d)-0.5)*2;
% mu=(dom(:,1)'+dom(:,2)')/2;
% h=-dom(:,1)'+dom(:,2)';
% for i=1:1:d
%     Xrnd(:,i)=(h(i)/2)*Xrnd(:,i)+mu(i);
% end



[Xrnd,wrnd]=smolyak_sparse_grid_modf(zeros(d,1),eye(d),d,4,'GLgn');
Nrnd=length(wrnd);
mu=(dom(:,1)'+dom(:,2)')/2;
h=-dom(:,1)'+dom(:,2)';
for i=1:1:d
    Xrnd(:,i)=(h(i)/2)*Xrnd(:,i)+mu(i);
end

% [Xrnd2,wrnd2]=smolyak_sparse_grid_modf(zeros(d,1),eye(d),d,4,'GLgn');
% Nrnd2=length(wrnd2);
% mu=(dom(:,1)'/3+dom(:,2)'/3)/2;
% h=-dom(:,1)'/3+dom(:,2)'/3;
% for i=1:1:d
%     Xrnd2(:,i)=(h(i)/2)*Xrnd2(:,i)+mu(i);
% end
% Xrnd=[Xrnd;Xrnd2];



Nrnd=size(Xrnd,1);

% % the lqr initial cost at these  collocation points
S=eye(d);
V0=zeros(length(w),1);
for i=1:1:length(w)
    V0(i)=X(i,:)*S*X(i,:)';
end


%% generating all the polynomia basis and the correponding matrices here

order_poly_approx=8;

if exist('PHI_6D_8M.m','file')==2
    phi=@(x)PHI_6D_8M(x);
    Jac=@(x)JAC_6D_8M(x);
else
    
    Pbasis=Basis_polyND(d,order_poly_approx);
    %     Pbasis=RemBasis_polyND(Pbasis,[1,3,5,7]);
    
    ss=GenMfile_MatrixOfPolys(Pbasis,'PHI_6D_8M','');
    phi = str2func(strcat('@(x)',ss));
    
    Jac=Jacobian_VectorOfPolys(Pbasis);
    ss=GenMfile_MatrixOfPolys(Jac,'JAC_6D_8M','');
    Jac = str2func(strcat('@(x)',ss));
    
end





m=length(phi(X(1,:)));

% PPT=AddSubMultiply_MatrixOfPolys(Pbasis,Pbasis','multiply');
% ss=GenMfile_MatrixOfPolys(PPT,'','');
% PPT = str2func(strcat('@(x)',ss));


Xdd=[-1.6000   -1.6000    1.6000    0.8000         0         0
   -1.6000    1.6000   -1.6000    0.8000         0         0
    1.6000   -1.6000   -1.6000   -0.8000         0         0
    1.6000    1.6000    1.6000   -0.8000         0         0
   -1.6000   -1.6000   -1.6000    0.8000    0.8000         0
   -1.6000   -1.6000    1.6000    0.8000    0.8000         0
   -1.6000    1.6000   -1.6000    0.8000   -0.8000         0
   -1.6000    1.6000    1.6000    0.8000   -0.8000         0
    1.6000   -1.6000   -1.6000   -0.8000    0.8000         0
    1.6000   -1.6000    1.6000   -0.8000    0.8000         0
    1.6000    1.6000   -1.6000   -0.8000   -0.8000         0
    1.6000    1.6000    1.6000   -0.8000   -0.8000         0];


XA=[];

XA=[general_conj_axis(d,d);general_conj_axis(d,2)];
for i=1:1:d
    XA(XA(:,i)==-1,i)=dom(i,1);
    XA(XA(:,i)==1,i)=dom(i,2);
end
XA=[X;XA;Xdd];

A=zeros(size(XA,1),m);
for i=1:1:size(XA,1)
    A(i,:)=phi(XA(i,:));
end

Abnd=zeros(N,m);
for i=1:1:N
    A(i,:)=phi(X(i,:));
end



BmapType='interp';
% BmapType='2norm';

if strcmpi(BmapType,'2norm')
    % check for simple quadratic
    [XX,ww] = GLeg_pts(7*ones(1,d),bdd_low,bdd_up);
    M=0;
    for i=1:1:length(ww)
        M=M+ww(i)*phi(XX(i,:))*phi(XX(i,:))';
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



type='CcoeffNOiter';
% type='Ccoeff';
% type='Vcoeff';
% type='Galerkin';

Csave=[];

for iter=1:1:2
    
    
    
    
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
            F=vertcat(F,     -xi'*Q*xi-u(xi)'*R*u(xi)    ) ;
        end
        
        %Vc=D\(F);
        %  c=Apinv*Vc;
        
        % derivarive constraint
        %                 D=vertcat(D,Jac([0,0])');
        %                 F=vertcat(F,zeros(d,1));
        
        %     % 0 at the origin constraint
        D=vertcat(D,phi(zeros(1,d))');
        F=vertcat(F,0);
        
        %adding the inequality constraints
        
        %                 keyboard
        
        Deq=D;
        Feq=F;
        indc=[];
        Np=5;
        TT=linspace(1,m-N,Np);
        TT=round(TT);
        TT=[1,TT];
        for jj=2:1:length(TT)
            
            for jk=1:1:(TT(jj)-TT(jj-1))
                bb=c;
                bb(indc)=1e15;
                [~,ic]=min(abs(bb));
                indc=unique(horzcat(indc,ic));
            end
            Z=zeros(length(indc),m);
            for hh=1:1:length(indc)
                Z(hh,indc(hh))=1;
            end
            Deq=vertcat(D,Z);
            Feq=vertcat(F,zeros(length(indc),1));
            
            indc
            
            cvx_begin
            variable c(m)
            minimize( norm( c,1 ) )
            subject to
            Deq * c == Feq
            A*c>=0
            cvx_end
            
            %             keyboard
            
        end
        c(indc)=0;
        %         keyboard
        
        %         options = optimset('Display','iter','Jacobian','off','Algorithm','interior-point','MaxFunEvals',1e7,'MaxIter',1e7);
        %         c=fmincon(@(c)sum(c.^2),c,Aineq,zeros(size(Aineq,1),1),D,F,[],[],[],options)
        %
    elseif strcmpi(type,'CcoeffNOiter')
        
        D=zeros(N,m);
        F=zeros(N,1);
        Aineq=zeros(Nrnd,m);
        for i=1:1:N
            xi=X(i,:)';
            Jxi=Jac(X(i,:));
            
            D(i,:)=     (f(xi)+g(xi)*u(xi))'*Jxi'   ;
            F(i)=-xi'*Q*xi-u(xi)'*R*u(xi)     ;
        end
        for i=1:1:Nrnd
            xrnd=Xrnd(i,:)';
            Aineq(i,:)=(f(xrnd)+g(xrnd)*u(xrnd))'*Jac(xrnd)';
        end
        
        Deq=D;
        Feq=F;
        
        %         c = quadprog(eye(m),zeros(m,1),Aineq,zeros(size(Aineq,1),1),Deq,Feq,[],[])
        %         c=linprog(ones(m,1),Aineq,zeros(size(Aineq,1),1),Deq,Feq,[],[]);
        cvx_begin
        variable c(m)
        minimize( norm( c,1 ) )
        subject to
        Deq * c == Feq
        A*c>=0
        %         Aineq*c<=0
        cvx_end
        
        %                     keyboard
        
        
        
    elseif strcmpi(type,'Galerkin')
        keyboard
        
        [Xint,wint] = GLeg_pts(8*ones(1,d),dom(:,1)',dom(:,2)');
        tic
        D=zeros(m);
        F=zeros(m,1);
        
        for j=1:1:length(wint)
            xi=Xint(j,:)';
            Jxi=Jac(xi');
            
            D=D+  wint(j)*repmat( (f(xi)+g(xi)*u(xi))'*Jxi',m,1).*repmat(phi(xi),1,m) ;
            F=F+  wint(j)*( -xi'*Q*xi-u(xi)'*R*u(xi) )*phi(xi)    ;
            j
        end
        D=D*prod(dom(:,2)-dom(:,1));
        F=F*prod(dom(:,2)-dom(:,1));
        
        D=[D;[1,zeros(1,m-1)]];
        F=[F;0];
        toc
        
        
        keyboard
        c=D\F;
        %         c=pinv(D)*F;
        
    end
    
    
    u=@(x)-0.5*inv(R)*g(x)'*Jac(x)'*c;
    Csave=vertcat(Csave,c');
    
    %     figure
    %     [xx,yy]=meshgrid(linspace(dom(1,1),dom(1,2),10)*0.8,linspace(dom(2,1),dom(2,2),10)*0.8);
    %     for i=1:1:size(xx,1)
    %         for j=1:1:size(xx,2)
    %             disp([i,j])
    %             [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);dom(3,2)*0.8;dom(4,2)*0.8;dom(5,2)*0.8;dom(6,2)*0.8]);
    %
    %             plot(x(:,1),x(:,2),'linewidth',2)
    %             hold on
    %             pause(0.5)
    %         end
    %     end
    %     xlabel('p_1')
    %     xlabel('p_2')
    %     plot_prop_paper
    %     if saveonoff==1
    %         saveas(gca,strcat(figname,'IterTrajs_',num2str(iter)),'eps')
    %         saveas(gca,strcat(figname,'IterTrajs_',num2str(iter)),'png')
    %         saveas(gca,strcat(figname,'IterTrajs_',num2str(iter)),'fig')
    %     end
    %
    
    
end
%%
ss=num2str(round(clock()));
ss(ss==' ')='';


Ngrid=4;
% Xgrid1=linspace(dom(1,1),dom(1,2),Ngrid);
% Xgrid2=linspace(dom(2,1),dom(2,2),Ngrid);
% Xgrid3=linspace(dom(3,1),dom(3,2),Ngrid);

[X1G,X2G,X3G,X4G,X5G,X6G]=ndgrid(linspace(dom(1,1),dom(1,2),Ngrid)*0.8  ,linspace(dom(2,1),dom(2,2),Ngrid)*0.8  ,linspace(dom(3,1),dom(3,2),Ngrid)*0.8,linspace(dom(4,1),dom(4,2),Ngrid)*0.8,linspace(dom(5,1),dom(5,2),Ngrid)*0.8,linspace(dom(6,1),dom(6,2),Ngrid)*0.8);
V0surf=zeros(size(X1G));
VFsurf=zeros(size(X1G));
Errsurf=zeros(size(X1G));

errnorm=zeros(prod(size(X1G)),1);
Vfinal=zeros(prod(size(X1G)),1);
Vinitial=zeros(prod(size(X1G)),1);
Xg=zeros(prod(size(X1G)),d);

kk=1;
NNN=size(X1G,1);
parfor i1=1:NNN
    HJB_6D_AttCtleg_PARFOR(i1,X1G,X2G,X3G,X4G,X5G,X6G,f,g,u0,Q,R,V0surf,VFsurf,Errsurf,d,c,phi,Jac)
end

for i1=1:NNN
    HH=load(strcat('HJBAttCtrl_',num2str(i1)));
    V0surf=V0surf+HH.V0surf;
    VFsurf=VFsurf+HH.VFsurf;
    Errsurf=Errsurf+HH.Errsurf;
end


kk=1;
for i1=1:1:size(X1G,1)
    for i2=1:1:size(X1G,2)
        for i3=1:1:size(X1G,3)
            for i4=1:1:size(X1G,4)
                for i5=1:1:size(X1G,5)
                    for i6=1:1:size(X1G,6)
                        
                        x1=X1G(i1,i2,i3,i4,i5,i6);
                        x2=X2G(i1,i2,i3,i4,i5,i6);
                        x3=X3G(i1,i2,i3,i4,i5,i6);
                        x4=X4G(i1,i2,i3,i4,i5,i6);
                        x5=X5G(i1,i2,i3,i4,i5,i6);
                        x6=X6G(i1,i2,i3,i4,i5,i6);
                        
                        xi=[x1;x2;x3;x4;x5;x6];
                        
                        Xg(kk,:)=xi;
                        
                        Vinitial(kk)=V0surf(i1,i2,i3,i4,i5,i6);
                        errnorm(kk)=Errsurf(i1,i2,i3,i4,i5,i6);
                        Vfinal(kk)=VFsurf(i1,i2,i3,i4,i5,i6);
                        
                        kk=kk+1
                    end
                end
            end
            
        end
    end
end
HJBerr=sqrt(mean(errnorm));
Xnorm=sqrt(sum(Xg.^2,2));
[Xnorm,indsort]=sort(Xnorm);
Vinitial=Vinitial(indsort);
Vfinal=Vfinal(indsort);

XXG=Xg(indsort,:);
XXG(Vfinal<0,:)


figure
plot(Xnorm,Vinitial,'ro-','linewidth',2)
hold on
plot(Xnorm,Vfinal,'bo-','linewidth',2)
legend('Initial','final')
xlabel('||x||')
ylabel('V(x)')
plot_prop_paper
if saveonoff==1
    saveas(gca,strcat(figname,'VFV0Xnorm_'),'eps')
    saveas(gca,strcat(figname,'VFV0Xnorm_'),'png')
    saveas(gca,strcat(figname,'VFV0Xnorm_'),'fig')
end

% for p=1:1:size(XG,3)
%     figure
%     mesh(XG(:,:,p),YG(:,:,p),V0surf(:,:,p))
%     hold on
%
%     surf(XG(:,:,p),YG(:,:,p),VFsurf(:,:,p))
%     alpha 0.5
%     xlabel('w_1')
%     ylabel('w_2')
%     plot_prop_paper
%     if saveonoff==1
%         saveas(gca,strcat(figname,'VFV0surf_',num2str(p)),'eps')
%         saveas(gca,strcat(figname,'VFV0surf_',num2str(p)),'png')
%     end
%
%     figure
%     contour(XG(:,:,p),YG(:,:,p),VFsurf(:,:,p),30)
%     Fx=YG.*ZG*(Iy-Iz)/Ix;
%     Fy=XG.*ZG*(Iz-Ix)/Iy;
%     Fz=XG.*YG*(Ix-Iy)/Iz;
%
%
%     hold on
%     quiver(XG(:,:,p),YG(:,:,p),Fx(:,:,p),Fy(:,:,p),1.5,'linewidth',2)
%     alpha 0.5
%     xlabel('w_1')
%     ylabel('w_2')
%     colorbar
%     plot_prop_paper
%     if saveonoff==1
%         saveas(gca,strcat(figname,'VFcont_',num2str(p)),'eps')
%         saveas(gca,strcat(figname,'VFcont_',num2str(p)),'png')
%     end
%
% end

%%

figure
[xx,yy]=meshgrid(linspace(dom(1,1),dom(1,2),10)*0.8,linspace(dom(2,1),dom(2,2),10)*0.8);
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        disp([i,j])
%         if i==10
%             plot(xx(i,j),yy(i,j),'ro','linewidth',2)
%         else
%             continue
%         end
        [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);dom(3,2)*0.8;dom(4,2)*0.8;dom(5,2)*0.8;dom(6,2)*0.8]);
        if norm(max(abs(x),[],1))<10
            plot(x(:,1),x(:,2),'linewidth',2)
        else
            plot(xx(i,j),yy(i,j),'ro','linewidth',2)
        end
        hold on
        pause(0.5)
    end
end
xlabel('p_1')
ylabel('p_2')
plot_prop_paper
if saveonoff==1
    saveas(gca,strcat(figname,'FinalTrajs'),'eps')
    saveas(gca,strcat(figname,'FinalTrajs'),'png')
    saveas(gca,strcat(figname,'FinalTrajs'),'fig')
end




save(strcat('HJBcollSol6D_',ss))




