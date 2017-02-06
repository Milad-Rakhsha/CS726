function [inform,xnew] = LBFGS(fun,x,lbfgsparams)
global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;
x.g = feval(fun,x.p,2);
x.f = feval(fun,x.p,1);
alfa=1.0;
params=struct('c1',1e-4,'c2',0.4,'maxit',20);
%%% Algorithm 7.5 - Numerical Optimization book (Second Edition)
I_Mat=eye(length(x.p));
xnew=x;
m=lbfgsparams.m;
%Most recent values (k-1) will be stored at the rightmost column
yk_Mat=zeros(length(x.p),m);
sk_Mat=zeros(length(x.p),m);
yk=zeros(length(x.p),1);
sk=zeros(length(x.p),1);
rho_vec=zeros(1,m);
alpha_i=zeros(1,m);
for k=1:lbfgsparams.maxit
    x=xnew;
    gradf=x.g;
    gamma=(yk'*sk)/(yk'*yk);
    if (~isfinite(gamma))
        gamma=0.001;
    end
    H0k=gamma*I_Mat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calc pk=-Hk*gradf from Algorithm 7.4
    for i=1:m % i =k-1, k-2, ... , k-m
        alpha_i(m-i+1)=rho_vec(m-i+1)*sk_Mat(:,m-i+1)'*gradf;
        gradf=gradf-alpha_i(m-i+1)*yk_Mat(:,m-i+1);
    end
    r=H0k*gradf;
    for i=1:m % i =k-m, k-m+1, ... , k-1
        beta=rho_vec(i)*yk_Mat(:,i)'*r;
        r=r+sk_Mat(:,i)*(alpha_i(i)-beta);
    end
    %stop with result Hk*gradfk=r
    pk=-r;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Solve a line search for alphak
    [alpha, xnew] = StepSize(fun, x, pk, alfa, params);
    %%% Having xnew and gradfnew in hand calc yk and sk
    yk=xnew.g-x.g;
    sk=xnew.p-x.p;
    %Discard the vector pair {s(k-m), y(k-m)} from storage
    yk_Mat(:,1)=[];
    sk_Mat(:,1)=[];
    % Append the new values to the end
    yk_Mat=[yk_Mat yk];
    sk_Mat=[sk_Mat sk];
    if(isfinite(1/(yk'*sk)))
        rho_vec(1)=[];
        rho_vec=[rho_vec  1/(yk'*sk)];
    else
        % This else should not be executed but just to make sure
        rho_vec(1)=[];
        rho_vec=[rho_vec  0];
    end
    %%%  Calculate Residual
    iter=iter+1;
    res=norm(xnew.g,inf) - lbfgsparams.toler*(1+abs(xnew.f));
    if (res<0)
        status=1;
        break
    end
end
xnew = struct('p',xnew.p,'f',xnew.f,'g',xnew.g);
inform = struct('status',status,'iter',k);
