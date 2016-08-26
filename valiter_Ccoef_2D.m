% value iteration by collocation
clear
clc
clear
close all

clear
%% Given system is
Q=diag([50,50]);
R=diag([5]);

figname='2DVanderpoll_lowCtrlCost_';
saveonoff=1;


d=2;

ep=1;
g=@(x)[0;1];
f=@(x)[x(2);-x(1)+x(2)-ep*x(1)^2*x(2)];


dom=[-3,3;
    -3,3];


%%

% -1.2 and -5
K=[-1.2,-1.2];
u=@(x)K*x(:);
u0=u;


% [xx,yy]=meshgrid(linspace(dom(1,1),dom(1,2),10));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j)]);
%         
%         plot(x(:,1),x(:,2),'linewidth',2)
%         hold on
%         pause(0.1)
%     end
% end
% xlabel('x_1')
% ylabel('x_2')
% plot_prop_paper
% if saveonoff==1
%     saveas(gca,strcat(figname,'InitialTrajs'),'eps')
%     saveas(gca,strcat(figname,'InitialTrajs'),'png')
% end



%% the points that will be used for collocation


Ncol_gh=4;

% bdd_low=-dom*ones(1,d);
% bdd_up=dom*ones(1,d);
X=[];
w=[];
X2=[];
w2=[];
% [X,w]=uniform_sigma_pts(dom(:,1)',dom(:,2)',6);
% [X2,w2]=uniform_sigma_pts(dom(:,1)'/3,dom(:,2)'/3,4);
[X,w] = GLeg_pts(Ncol_gh*ones(1,d),dom(:,1)',dom(:,2)');
[X2,w2] = GLeg_pts((Ncol_gh-1)*ones(1,d),dom(:,1)'/3,dom(:,2)'/3);

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



type='Ccoeff';
% type='Vcoeff';
% type='Galerkin';

Csave=[];

for iter=1:1:40
    
    
    
    
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
            
            
            
        end
        c(indc)=0;
        %         keyboard
        
        %         options = optimset('Display','iter','Jacobian','off','Algorithm','interior-point','MaxFunEvals',1e7,'MaxIter',1e7);
        %         c=fmincon(@(c)sum(c.^2),c,Aineq,zeros(size(Aineq,1),1),D,F,[],[],[],options)
        %
        
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
    
    
    u=@(x)-0.5*inv(R)*g(x)'*Jac(x)'*c;
    Csave=vertcat(Csave,c');
    
    
    
    
end
%%
ss=num2str(round(clock()));
ss(ss==' ')='';







Ngrid=20;
[Xgridx,Xgridy]=meshgrid(linspace(dom(1,1),dom(1,2),Ngrid)*0.8,linspace(dom(2,1),dom(2,2),Ngrid)*0.8);

V0surf=zeros(size(Xgridx));
VFsurf=zeros(size(Xgridx));
Errsurf=zeros(size(Xgridx));

errnorm=zeros(size(Xgridx,1),1);
Vfinal=zeros(size(Xgridx,1),1);
Vinitial=zeros(size(Xgridx,1),1);
Xg=zeros(size(Xgridx,1),d);
kk=1;
for i=1:1:size(Xgridx,1)
    for j=1:1:size(Xgridx,2)
            
            xi=[Xgridx(i,j);Xgridy(i,j)];
            Xg(kk,:)=xi;
            Vfinal(kk)=c'*phi(xi);
            [t,x]=ode45(@(t,x)[f(x(1:d))+g(x(1:d))*u0(x(1:d));x(1:d)'*Q*x(1:d)+u0(x(1:d))'*R*u0(x(1:d))],[0,200],[xi;0]);
            Vinitial(kk)=x(end,d+1);
            
            
            errnorm(kk)=HJBerror(Q,R,f,g,@(x)Jac(x),xi,c)^2;
            
            
            V0surf(i,j)=Vinitial(kk);
            VFsurf(i,j)=Vfinal(kk);
            Errsurf(i,j)=errnorm(kk);
            kk=kk+1;
    end
end
HJBerr=sqrt(mean(errnorm))

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
    saveas(gca,strcat(figname,'VFV0Xnorm'),'eps')
    saveas(gca,strcat(figname,'VFV0Xnorm'),'png')
end



figure

mesh(Xgridx,Xgridy,V0surf)
hold on
surf(Xgridx,Xgridy,VFsurf)
alpha 0.5
xlabel('x_1')
ylabel('x_2')
plot_prop_paper
if saveonoff==1
     saveas(gca,strcat(figname,'VFV0surf'),'eps')
        saveas(gca,strcat(figname,'VFV0surf'),'png')

end

figure
contour(Xgridx,Xgridy,VFsurf,linspace(0.01,400,20),'linewidth',2)
Fx=Xgridy;
Fy=-Xgridx+Xgridy-ep*Xgridx.^2.*Xgridy;
hold on
quiver(Xgridx,Xgridy,Fx,Fy,2,'linewidth',2)
alpha 0.5
colorbar
xlabel('x_1')
ylabel('x_2')
plot_prop_paper
if saveonoff==1
    saveas(gca,strcat(figname,'VFcont'),'eps')
    saveas(gca,strcat(figname,'VFcont'),'png')
end


figure
[xx,yy]=meshgrid(linspace(dom(1,1)*0.8,dom(1,2)*0.8,10));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        disp([i,j])
        [t,x]=ode45(@(t,x)f(x)+g(x)*u(x),[0,200],[xx(i,j);yy(i,j)]);
        
        plot(x(:,1),x(:,2),'linewidth',2)
        hold on
        pause(0.1)
    end
end
xlabel('x_1')
ylabel('x_2')
plot_prop_paper
if saveonoff==1
    saveas(gca,strcat(figname,'FinalTrajs'),'eps')
    saveas(gca,strcat(figname,'FinalTrajs'),'png')
end



save(strcat('HJBcollSol_2D_vanderpoll_',ss))