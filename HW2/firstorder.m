clear
clc
mu=0.01; L=1; kappa=L/mu;
n=100;
A = randn(n,n); [Q,R]=qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D = mu + D*(L-mu);
A = Q'*diag(D)*Q;


RES=1e-6;
iter_sd=0;
iter_sde=0;
iter_nest=0;
iter_conj=0;

for ITER=1:10
    x0 = randn(n,1); % use a different x0 for each trial
    
    if ITER==1
        res_1=[];
        res_2=[];
        res_3=[];
        res_4=[];
        
    end
    %    Steepest descent with k1=L
    res=1;
    x=x0;
    while (res>RES)
        x_old=x;
        gradf=A*x_old;
        alpha=1/L;
        x=x_old-gradf*alpha;
        res=abs(1/2*transpose(x)*A*x);
        iter_sd=iter_sd+1;
        if ITER==1
            res_1(iter_sd,:)=[iter_sd,res];
        end
    end
    %     Steepest descent with exact line search.
    res=1;
    x=x0;
    while (res>RES)
        x_old=x;
        gradf=A*x_old;
        alpha=1*(transpose(x)*A*gradf)/(transpose(gradf)*A*gradf);
        x=x_old-gradf*alpha;
        res=abs(1/2*transpose(x)*A*x);
        iter_sde=iter_sde+1;
        if ITER==1
            res_2(iter_sde,:)=[iter_sde,res];
        end
    end
    
    %     Nesterov's optimal method
    res=1;
    x=x0;
    x_old=x0;
    x_inter=x0;
    while (res>RES)
        x_old=x_inter;
        beta=(sqrt(kappa)-1)/(sqrt(kappa)+1);
        alpha=1/L;
        yk=x+beta*(x-x_old);
        gradf=A*(yk);
        x_inter=x;
        x=yk-alpha*gradf;
        res=abs(1/2*transpose(x)*A*x);
        iter_nest=iter_nest+1;
        if ITER==1
            res_3(iter_nest,:)=[iter_nest,res];
        end
    end
    
    
    %     Conjugate gradient method from p.108 of the text.
    res=1;
    x=x0;
    y=x0;
    r=A*x0;
    p=-r;
    beta=(r'*A*y)/(y'*A*y);
    while (res>RES)
        x_old=x;
        alpha=-r'*p/(p'*A'*p);
        x=x_old+alpha*p;
        r=A*x;
        beta=(r'*A*p)/(p'*A*p);
        p=-r+beta*p;
        res=abs(1/2*transpose(x)*A*x);
        iter_conj=iter_conj+1;
        if ITER==1
            res_4(iter_conj,:)=[iter_conj,res];
        end
        
    end
    
end


% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1000, 600]);
% plot((res_1(:,1)),log10(res_1(:,2)),'k','lineWidth',2);
% hold on
% plot((res_2(:,1)),log10(res_2(:,2)),'b--','lineWidth',2);
% plot((res_3(:,1)),log10(res_3(:,2)),'r.-.','lineWidth',2);
% plot((res_4(:,1)),log10(res_4(:,2)),'m-o','lineWidth',2);
% grid on
% 
% xlabel('Iteration','FontSize',18,'fontName','Times New Roman');
% ylabel('log (f(x)-f(x^{*}))','FontSize',18,'fontName','Times New Roman')
% set(gca,'FontSize',18)
% AX=legend('SD:const','SD:exact','Nesterov','CG');
% set(AX,'FontSize',20,'fontName','Times New Roman');
% set(FigHandle,'Units','Inches');
% pos = get(FigHandle,'Position');
% set(FigHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(FigHandle,'Iter_Res.pdf','-dpdf','-r0')


fprintf(1,' steepest descent - fixed steps : %7.1f\n', iter_sd/10);
fprintf(1,' steepest descent - exact steps : %7.1f\n', iter_sde/10);
fprintf(1,' Nesterov : %7.1f\n', iter_nest/10);
fprintf(1,' conjugate gradient : %7.1f\n', iter_conj/10);
