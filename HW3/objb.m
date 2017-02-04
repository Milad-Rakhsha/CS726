function varargout = objb(x,mode)

global numf numg numH

argout = 0;
if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  varargout(argout) = {x(1)^2 + 5*x(1)*x(2) + 100*x(2)^2 - x(1) + 4*x(2)};
end
if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  varargout(argout) = {[2*x(1) + 5*x(2) - 1; 5*x(1) + 200*x(2) + 4]};
end
if bitand(mode,4) 
  numH = numH + 1;
  argout = argout + 1;
  varargout(argout) = {[2 5; 5 200]};
end
return;
