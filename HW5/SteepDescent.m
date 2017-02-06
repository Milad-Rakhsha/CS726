function [inform,x] = SteepDescent(fun,x,sdparams)

global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;
% x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
x.f = feval(fun,x.p,1);
pk=-x.g;
alfa=1;
params=struct('c1',0.01,'c2',0.85,'maxit',100);

for k=1:sdparams.maxit
    [alfa, xnew] = StepSize(fun, x, pk, alfa, params);
    x=xnew;
    iter=iter+1;
    pk=-x.g;
    res=norm(xnew.g,inf) - sdparams.toler*(1+abs(xnew.f));
    if (res<0)
        status=1;
        break
    end
end

x.g=feval(fun,x.p,2);
x = struct('p',x.p,'f',x.f,'g',x.g);
inform = struct('status',status,'iter',k);