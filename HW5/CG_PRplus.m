function [inform,x] = CG_PRplus(fun,x,sdparams)


global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;
% x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
x.f = feval(fun,x.p,1);
pk=-x.g;
alfa=1.0;
params=struct('c1',0.01,'c2',0.3,'maxit',100);
xnew=x;
for k=1:sdparams.maxit
    x=xnew;
    %%%Solve a line search
    gradf=x.g;
    [alpha, xnew] = StepSize(fun, x, pk, alfa, params);
    gradfnew=xnew.g;
    %%%Update
    beta=(gradfnew'*(gradfnew-gradf))/(gradf'*gradf);
    beta=max(0,beta);
    pk=-gradfnew+beta*pk;
    iter=iter+1;
    res=norm(xnew.g,inf) - sdparams.toler*(1+abs(xnew.f));
    if (res<0)
        status=1;
        break
    end
end

x.g=feval(fun,x.p,2);
x = struct('p',x.p,'f',x.f,'g',x.g);
inform = struct('status',status,'iter',k);