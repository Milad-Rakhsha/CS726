function varargout = tridia(x,mode)

  global numf  numg
  
argout = 0;
n=size(x,1);
if bitand(mode,1) 
  argout = argout + 1;
  f = (1.0-x(1))^2;
  for i=1:n-1
    f=f+(x(i)-2*x(i+1))^4;
  end
  f = 0.5*f;
  numf=numf+1;
  varargout(argout) = {f};
end
if bitand(mode,2) 
  argout = argout + 1;
  grad = zeros(n,1);
  grad(1)=(x(1)-1.0) + 2*(x(1)-2*x(2))^3;
  for i=2:n-1
    grad(i) = 2*(x(i)-2*x(i+1))^3 - 4*(x(i-1)-2*x(i))^3;
  end
  grad(n) = -4*(x(n-1)-2*x(n))^3;
  numg = numg+1;
  varargout(argout) = {grad};
end
return;
