function varargout = obja(x,mode)

global numf numg numH

argout = 0;
if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  varargout(argout) = {x(1)^2 + 5*x(2)^2 + x(1) - 5*x(2)};
end
if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  varargout(argout) = {[2*x(1) + 1; 10*x(2) - 5]};
end
if bitand(mode,4) 
  numH = numH + 1;
  argout = argout + 1;
  varargout(argout) = {[2 0; 0 10]};
end
return;
