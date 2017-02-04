function [inform,x] = DoglegTR(fun, x, trparams)

global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;
x.f = feval(fun,x.p,1);

% Algorithm 4.1 parameters
dhat =trparams.hatDelta; % must be greater than 0
dk = trparams.delta; % initial delta, between 0 and dhat
eta = trparams.eta; % [0,0.25)

% Perform interation
for k=1:trparams.maxit
    % Choose pk using Newton's direction
    gk = feval(fun,x.p,2); %gk=gk/norm(gk,2);  % normalize gk
    Bk = feval(fun,x.p,4);
    pk = -dk*(Bk\gk); % multiply by dk so that the length of pk = dk;
    
    % Check whether the hessian is approximately positive definite at xk
    % Actually it is supposed to be p'*B*p >0 for all p in the trust region
    if pk.'*Bk*pk <= 0
        error('Hessian is not positive definite at xk')
    end
    
    % Form rhok
    fk = feval(fun,x.p,1); fkpk =feval(fun,x.p+pk,1);
    mk0 = fk; mkp = fk + gk.'*pk + 0.5*pk.'*Bk*pk;
    rhok = (fk - fkpk)/(mk0 - mkp);
    
    if rhok < 0.25
        dk = 0.25*dk;
    else
        if rhok > 0.75 && abs(norm(pk)-dk) < eps
            dk = min(2*dk,dhat);
        end
    end
    
    if rhok > eta
        % update xk to xk+1
        xold=x.p;
        x.p = x.p + pk;
    end
            iter=iter+1;
        res=norm(feval(fun,x.p,2));
        if res < trparams.toler
            status=1;
            break;
        end
    
end

x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
x.h = feval(fun,x.p,4);

x = struct('p',x.p,'f',x.f,'g',x.g,'h',x.h);
inform = struct('status',status,'iter',iter);



