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
syms xs;
% Perform interation
for k=1:trparams.maxit
    
    g = feval(fun,x.p,2);
    B = feval(fun,x.p,4);
    
    if g'*B*g <= 0
        [Evec, EV_matrix]=eig(B);
        EV=diag(EV_matrix);
        lambda_i_new=max(trparams.delta*ones(length(x.p),1),EV);
        B=Evec*(EV_matrix+diag(lambda_i_new))*Evec';
        eig(B);
        taw=1;
    else
        g_norm=norm(g,2);
        taw=min(g_norm^3/(dk*g'*B*g),1);
    end
    
    pu_taw=taw*-g/g_norm*dk;
    pu = -g'*g/(g'*B*g)*g;
    pb=-inv(B)*g;
    
    if (norm(pb,2)<dk)
        p=pb;
        
    else
        %have to solve for taw=x+1
        % ||pu+x*(pb-pu)||=dk
        pu=pu_taw;
        a=(pb-pu)'*(pb-pu);
        b=2*pu'*(pb-pu);
        c=pu'*pu-dk^2;
        xsol=(-b+sqrt(b^2-4*a*c))/2/a;
        taw=xsol+1;
        p=pu+(taw-1)*(pb-pu);
    end
    
    % Form rhok
    fk = feval(fun,x.p,1); fkpk =feval(fun,x.p+p,1);
    mk0 = fk; mkp = fk + g'*p + 0.5*p'*B*p;
    rhok = (fk - fkpk)/(mk0 - mkp);
    
    if rhok < 0.25
        dk = 0.25*dk;
    else
        if rhok > 0.75 && abs(norm(p)-dk) < eps
            dk = min(2*dk,dhat);
        end
    end
    
    if rhok > eta
        x.p = x.p + p;
        iter=iter+1;
        res=norm(feval(fun,x.p,2));
    end
    
    if res < trparams.toler
        status=1;
        break
    end
    
end

x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
x.h = feval(fun,x.p,4);

x = struct('p',x.p,'f',x.f,'g',x.g,'h',x.h);
inform = struct('status',status,'iter',iter);



