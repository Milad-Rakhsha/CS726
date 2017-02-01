clear
clc
syms x1 x2 f
f=(x1+x2^2)^2;
gradf=[diff(f,x1); diff(f,x2)]

x=[1;0];
p=[-1;0];
res=1;
alpha=0.3;
while res>eps
    
    x_old=x;
    x=x_old-alpha*vpa(subs(gradf,{x1,x2},{x_old(1), x_old(2)}),8);
    res=vpa(norm(x_old-x,2),6);
    fval=vpa(subs(f,{x1,x2},{x(1), x(2)}),3);
    display(['residual is:', num2str(double(res)), ', f value is:',num2str(double(fval))])
    
end

vpa(subs(f,{x1,x2},{x(1), x(2)}))
x