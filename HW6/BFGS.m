function [inform,xnew] = BFGS(fun,x,qnparams)
global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;
x.g = feval(fun,x.p,2);
x.f = feval(fun,x.p,1);
alfa=1.0;
params=struct('c1',1e-4,'c2',0.4,'maxit',20);
%%% Algorithm 6.1 - Numerical Optimization book (Second Edition)
Hk=eye(length(x.p));
I_Mat=eye(length(x.p));
xnew=x;
Hnew=Hk;
for k=1:qnparams.maxit
    x=xnew;
    Hk=Hnew;
    gradf=x.g;
    %%% Calc pk from the previous Hk the approximate Inv Hessian
    pk=-Hk*gradf;
    %%% Solve a line search for alphak
    [alpha, xnew] = StepSize(fun, x, pk, alfa, params);
    %%% Having xnew and gradfnew in hand calc yk and sk
    yk=xnew.g-x.g;
    sk=xnew.p-x.p;
    %%% Update: Compute H_(k+1) based on Eq (6.17)
    rhok=1/(yk'*sk);
    
    if (iter==0)
        %%% We can do a better starting H0 than Identity, using Eq (6.20)
        Hk=(yk'*sk)/(yk'*yk)*I_Mat;
    end
    Hnew= (I_Mat-rhok*sk*yk')*Hk*(I_Mat-rhok*yk*sk')+ rhok*(sk)*sk';
    
    %%%  Calculate Residual
    iter=iter+1;
    res=norm(xnew.g,inf) - qnparams.toler*(1+abs(xnew.f));
    if (res<0)
        status=1;
        break
    end
end
xnew = struct('p',xnew.p,'f',xnew.f,'g',xnew.g);
inform = struct('status',status,'iter',k);