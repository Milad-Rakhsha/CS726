function varargout = nls_residb(x,mode)

  global numf numg
  
  argout = 0;
  a = [0.13294; -0.244378; 0.325895];
  d = [2.5074; -1.36401; 1.02282];
  H = [-0.564255 0.392417; -0.404979 0.927589; -0.0735084 0.535493];
  B = [5.66598 2.77141; 2.77141 2.12413];
  q = (x-ones(size(x)))'*B*(x-ones(size(x)));
  % d = d*100;
  d = d*50;
  r = a + H*x + 0.5*q*d;
  f = 0.5*sum(r.^2);
  if bitand(mode,1) 
    argout = argout + 1;
    numf = numf+1;
    varargout(argout) = {f};
  end
  if bitand(mode,2) 
    argout = argout + 1;
    J = H + d*(B*(x-ones(size(x))))';
    numg = numg+1;
    varargout(argout) = {J'*r};
  end
  return;
  