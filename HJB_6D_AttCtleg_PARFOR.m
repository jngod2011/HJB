function HJB_6D_AttCtleg_PARFOR(i1,X1G,X2G,X3G,X4G,X5G,X6G,f,g,u0,Q,R,V0surf,VFsurf,Errsurf,d,c,phi,Jac)

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
   

                        [t,x]=ode45(@(t,x)[f(x(1:d))+g(x(1:d))*u0(x(1:d));x(1:d)'*Q*x(1:d)+u0(x(1:d))'*R*u0(x(1:d))],[0,200],[xi;0]);
                        
                        V0surf(i1,i2,i3,i4,i5,i6)=x(end,d+1);
                        VFsurf(i1,i2,i3,i4,i5,i6)=c'*phi(xi);
                        Errsurf(i1,i2,i3,i4,i5,i6)=HJBerror(Q,R,f,g,@(x)Jac(x),xi,c)^2;
                        
                    end
                end
            end
            
        end
    end
    
    save(strcat('HJBAttCtrl_',num2str(i1)),'V0surf','VFsurf','Errsurf')
    
    
end