function err=HJBerror(Q,R,f,g,J,x,c)

err=x'*Q*x+c'*J(x)*f(x)-0.25*c'*J(x)*g(x)*inv(R)*g(x)'*J(x)'*c;