
clear
clc
global numf numg 

% use same parameters throughout
qnparams = struct('maxit',1000,'toler',1.0e-6);

% testing BFGS

x = struct('p',ones(3,1));

[inform,xnew] = BFGS(@nls_resida,x,qnparams);
fprintf('\n\n Function resida running BFGS\n');
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

x = struct('p',zeros(2,1));

[inform,xnew] = BFGS(@nls_residb,x,qnparams);
fprintf('\n\n Function residb running BFGS\n');
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


x = struct('p',ones(100,1));

[inform,xnew] = BFGS(@xpowsing,x,qnparams);
fprintf('\n\n Function xpowsing running BFGS\n');
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', qnparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
% fprintf('  Ending point: '); fprintf('%8.4g ',xnew.p);
fprintf('\n  Ending value: '); fprintf('%8.4g ',xnew.f);
fprintf('; No. function evaluations: %d',numf);
% fprintf('\n  Ending gradient: '); fprintf('%8.4g ',xnew.g);
fprintf('; No. gradient evaluations %d',numg);
fprintf('\n  Norm of ending gradient: %8.4g\n\n\n', norm(xnew.g));



n=1000;

mvec = [3 5 8 12 20]';

for i=1:5

  x = struct('p',ones(n,1));
  lbfgsparams=struct('toler',1.e-4,'maxit',1000,'m',mvec(i));

  [inform,xnew] = LBFGS(@tridia,x,lbfgsparams);
  fprintf('\n\n Function tridia running LBFGS with m=%d\n', lbfgsparams.m);
  if inform.status == 0
    fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
    fprintf('gradient size decreasing below %8.4g.\n', lbfgsparams.toler);
  else
    fprintf('Success: %d steps taken\n', inform.iter);
  end
  fprintf('\n  Ending value: '); fprintf('%8.4g ',xnew.f);
  fprintf('; No. function evaluations: %d',numf);
  fprintf('; No. gradient evaluations %d',numg);
  fprintf('\n  Norm of ending gradient: %8.4g\n\n\n', norm(xnew.g));
  
end

  
