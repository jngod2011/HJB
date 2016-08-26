function [f,J]=quad_equal_const(V,B,d,z)
V=V(:);
n=length(d);
f=zeros(n,1);
for i=1:1:n
    f(i)=V'*B{i}*V+V'*d{i}+z{i};
end

% computing the jacobian
J=zeros(n,length(V));
for i=1:1:n
    J(i,:)=(B{i}+B{i}')*V+d{i};
     
end















