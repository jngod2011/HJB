% value iteration by collocation
clear
clc
clear
close all
clear
%% Given system is
Q=diag([50,50,50]);
R=diag([3,3,3]);

figname='3DSpinStab_lowCtrlCost_fullcoll_';
saveonoff=0;

d=3;

Ix=14;
Iy=10;
Iz=8;

g=@(x)[diag([1/Ix,1/Iy,1/Iz])];
f=@(x)[x(2)*x(3)*(Iy-Iz)/Ix;
    x(1)*x(3)*(Iz-Ix)/Iy;
    x(1)*x(2)*(Ix-Iy)/Iz];


dom=[-1,1;
    -1,1;
    -1,1];



%%

% 1 and 5
K=[-1,0,0;
    0,-1,0;
    0,0,-1];
u=@(x)K*x(:);
u0=u;

[xx,yy]=meshgrid(linspace(dom(1,1),dom(1,2),10));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);1]);
% 
%         plot(x(:,1),x(:,2),'linewidth',2)
%         hold on
%         pause(0.1)
%     end
% end
% plot_prop_paper
% xlabel('w1')
% ylabel('w2')
% if saveonoff==1
%     saveas(gca,strcat(figname,'InitialTrajs'),'eps')
%     saveas(gca,strcat(figname,'InitialTrajs'),'png')
% end


%% the points that will be used for collocation


Ncol_gh=3;

% bdd_low=-dom*ones(1,d);
% bdd_up=dom*ones(1,d);
X=[];
w=[];
X2=[];
w2=[];
[X,w]=uniform_sigma_pts(dom(:,1)',dom(:,2)',6);
[X2,w2]=uniform_sigma_pts(dom(:,1)'/2,dom(:,2)'/2,4);
% [X,w] = GLeg_pts(Ncol_gh*ones(1,d),bdd_low,bdd_up);
% [X2,w2] = GLeg_pts((Ncol_gh-1)*ones(1,d),bdd_low/(2*dom),bdd_up/(2*dom));

X=[X;X2];
w=[w;w2];



% X=X+randn(size(X))
W=diag(w);

N=length(w);

% % the lqr initial cost at these  collocation points
S=eye(d);
V0=zeros(length(w),1);
for i=1:1:length(w)
    V0(i)=X(i,:)*S*X(i,:)';
end


%% generating all the polynomia basis and the correponding matrices here

order_poly_approx=6;

Pbasis=Basis_polyND(d,order_poly_approx);
ss=GenMfile_MatrixOfPolys(Pbasis,'','');
phi = str2func(strcat('@(x)',ss));


Jac=Jacobian_VectorOfPolys(Pbasis);
ss=GenMfile_MatrixOfPolys(Jac,'','');
Jac = str2func(strcat('@(x)',ss));

m=length(Pbasis);

% PPT=AddSubMultiply_MatrixOfPolys(Pbasis,Pbasis','multiply');
% ss=GenMfile_MatrixOfPolys(PPT,'','');
% PPT = str2func(strcat('@(x)',ss));

A=zeros(N,m);
for i=1:1:N
    A(i,:)=phi(X(i,:));
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



% type='Ccoeff';
% type='Vcoeff';
type='Galerkin';
% type='FullColl';
XX=[X;rand(m-N,d)*2-1]; 
Csave=[];

for iter=1:1:5
    iter
    
    
    
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
            minimize( norm( c, 1 ) )
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
        
    elseif strcmpi(type,'Galerkin')
        [Xint,wint] = GLeg_pts(8*ones(1,d),dom(:,1)',dom(:,2)');
        D=zeros(m);
        F=zeros(m,1);
        
        for j=1:1:length(wint)
            xi=Xint(j,:)';
            Jxi=Jac(xi');
            
            D=D+  wint(j)*repmat( (f(xi)+g(xi)*u(xi))'*Jxi',m,1).*repmat(phi(xi),1,m) ;
            F=F+  wint(j)*( -xi'*Q*xi-u(xi)'*R*u(xi) )*phi(xi)    ;
            
        end
        D=D*prod(dom(:,2)-dom(:,1));
        F=F*prod(dom(:,2)-dom(:,1));
        
        D=[D;[1,zeros(1,m-1)]];
        F=[F;0];
        
        
        
%     keyboard
          c=D\F;
%         c=pinv(D)*F;
        
       elseif strcmpi(type,'FullColl')
        
          
        D=[];
        F=[];
        Aineq=-A;
        for i=1:1:m
            xi=XX(i,:)';
            Jxi=Jac(XX(i,:));
            %             Aineq=vertcat(Aineq,f(xi)'*Jxi');
            D=vertcat(D,     (f(xi)+g(xi)*u(xi))'*Jxi'   ) ;
            F=vertcat(F,     -xi'*Q*xi-u(xi)'*R*u(xi)      ) ;
        end
        c=D\(F);
        
        
        
    end
    
    
    u=@(x)-0.5*inv(R)*g(x)'*Jac(x)'*c;
    Csave=vertcat(Csave,c');
    
    
    
    
end
%%

% keyboard


ss=num2str(round(clock()));
ss(ss==' ')='';





Ngrid=10;
% Xgrid1=linspace(dom(1,1),dom(1,2),Ngrid);
% Xgrid2=linspace(dom(2,1),dom(2,2),Ngrid);
% Xgrid3=linspace(dom(3,1),dom(3,2),Ngrid);

[XG,YG,ZG]=ndgrid(linspace(dom(1,1),dom(1,2),Ngrid)  ,linspace(dom(2,1),dom(2,2),Ngrid)  ,linspace(dom(3,1),dom(3,2),Ngrid));
V0surf=zeros(size(XG));
VFsurf=zeros(size(XG));
Errsurf=zeros(size(XG));

errnorm=zeros(prod(size(XG)),1);
Vfinal=zeros(prod(size(XG)),1);
Vinitial=zeros(prod(size(XG)),1);
Xg=zeros(prod(size(XG)),d);
kk=1;
for i1=1:1:size(XG,1)
    for i2=1:1:size(XG,2)
        for i3=1:1:size(XG,3)
            x1=XG(i1,i2,i3);
            x2=YG(i1,i2,i3);
            x3=ZG(i1,i2,i3);
            
            xi=[x1;x2;x3];
            Xg(kk,:)=xi;
            Vfinal(kk)=c'*phi(xi);
            [t,x]=ode45(@(t,x)[f(x(1:3))+g(x(1:3))*u0(x(1:3));x(1:3)'*Q*x(1:3)+u0(x(1:3))'*R*u0(x(1:3))],[0,200],[xi;0]);
            Vinitial(kk)=x(end,4);
            errnorm(kk)=HJBerror(Q,R,f,g,@(x)Jac(x),xi,c)^2;
            
            V0surf(i1,i2,i3)=Vinitial(kk);
            VFsurf(i1,i2,i3)=Vfinal(kk);
            Errsurf(i1,i2,i3)=errnorm(kk);
            
            kk=kk+1;
            
        end
    end
end
HJBerr=sqrt(mean(errnorm));

Xnorm=sqrt(sum(Xg.^2,2));
[Xnorm,indsort]=sort(Xnorm);
Vinitial=Vinitial(indsort);
Vfinal=Vfinal(indsort);


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
end

for p=1:1:size(XG,3)
    figure
    mesh(XG(:,:,p),YG(:,:,p),V0surf(:,:,p))
    hold on
    
    surf(XG(:,:,p),YG(:,:,p),VFsurf(:,:,p))
    alpha 0.5
    xlabel('w_1')
    ylabel('w_2')
    plot_prop_paper
    if saveonoff==1
        saveas(gca,strcat(figname,'VFV0surf_',num2str(p)),'eps')
        saveas(gca,strcat(figname,'VFV0surf_',num2str(p)),'png')
    end
    
    figure
    contour(XG(:,:,p),YG(:,:,p),VFsurf(:,:,p),30)
    Fx=YG.*ZG*(Iy-Iz)/Ix;
    Fy=XG.*ZG*(Iz-Ix)/Iy;
    Fz=XG.*YG*(Ix-Iy)/Iz;
    
    
    hold on
    quiver(XG(:,:,p),YG(:,:,p),Fx(:,:,p),Fy(:,:,p),1.5,'linewidth',2)
    alpha 0.5
    xlabel('w_1')
    ylabel('w_2')
    colorbar
    plot_prop_paper
    if saveonoff==1
        saveas(gca,strcat(figname,'VFcont_',num2str(p)),'eps')
        saveas(gca,strcat(figname,'VFcont_',num2str(p)),'png')
    end
    
end

%%
figure
[xx,yy]=meshgrid(linspace(dom(1,1)*0.9,dom(1,2)*0.9,10));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        disp([i,j])
        [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j);1.5]);
        
        plot(x(:,1),x(:,2),'linewidth',2)
        hold on
        pause(0.1)
    end
end
xlabel('w_1')
ylabel('w_2')
plot_prop_paper
if saveonoff==1
    saveas(gca,strcat(figname,'FinalTrajs'),'eps')
    saveas(gca,strcat(figname,'FinalTrajs'),'png')
end


save(strcat('HJBcollSol_3D_',ss))


