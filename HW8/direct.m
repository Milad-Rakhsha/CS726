function [inform,xnew] = direct(fun,x,directparams)
global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
iter=0;

% This is the dimension of the input
n=length(x.p);
% vector of search directions
D=ones(2*n,1);
% This is supposed to be {1,-1,1,-1,...} so,
D(2:2:end,1)=-1;
gamma=1;

%The following implements the Algorithm 9.2 from the book

for k=1:directparams.maxit
    
    if (gamma<directparams.toler)
        status=1;
        break
    end
    dir_IDX=randperm(2*n);
    x.f=feval(fun,x.p,1);
    

    for i=1:length(D)
        pk=zeros(n,1);
        ith=dir_IDX(i);
        pk(ceil(ith/2))=D(ith);
        %We are using \rho(\gamma)=\gamma^2
        if(feval(fun,x.p+gamma*pk,1)<x.f-gamma^2)
            xnew=x.p+gamma*pk;
            gammanew=directparams.phi*gamma;
            break;
        else
            xnew=x.p;
            gammanew=directparams.theta*gamma;
        end
        
    end
    
    x.p=xnew;
    gamma=gammanew;
    x.f=feval(fun,x.p,1);
    iter=iter+1;
    
end

xnew = struct('p',x.p,'f',x.f);
inform = struct('status',status,'iter',iter);