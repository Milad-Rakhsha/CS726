clear;clc
global numf numg 

% use same parameters throughout
qnparams = struct('maxit',1000,'toler',1.0e-6);

x0=randn(8,1);

x = struct('p',x0);

[inform,xnew] = BFGS(@xpowsing,x,qnparams);
fprintf('\n\n Function xpowsing running BFGS\n');
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', qnparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('  Ending point: '); fprintf('%8.4g ',xnew.p);
fprintf('\n  Ending value: '); fprintf('%8.4g ',xnew.f);
fprintf('; No. function evaluations: %d',numf);
% fprintf('\n  Ending gradient: '); fprintf('%8.4g ',xnew.g);
fprintf('; No. gradient evaluations %d',numg);
fprintf('\n  Norm of ending gradient: %8.4g\n\n\n', norm(xnew.g));


directparams = struct('maxit',100000,'toler',1.0e-6,'theta',.5,'phi',2);

x = struct('p',x0);
[inform,xnew] = direct(@xpowsing,x,directparams);
fprintf('\n\n Function xpowsing running Direct Search\n');
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', directparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('  Ending point: '); fprintf('%8.4g ',xnew.p);
fprintf('\n  Ending value: '); fprintf('%8.4g ',xnew.f);
fprintf('; No. function evaluations: %d\n\n\n',numf);
% fprintf('\n  Ending gradient: '); fprintf('%8.4g ',xnew.g);
%fprintf('; No. gradient evaluations %d',numg);
%fprintf('\n  Norm of ending gradient: %8.4g\n\n\n', norm(xnew.g));

