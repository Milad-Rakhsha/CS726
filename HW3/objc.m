function varargout = objc(x,mode)

global numf numg numH

argout = 0;
z= x(2) - x(1)^2;
if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  varargout(argout) = {100*z^2 + (1-x(1))^2};
end
if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  varargout(argout) = {[-400*x(1)*z + 2*x(1) - 2; 200*z]};
end
if bitand(mode,4) 
  numH = numH + 1;
  argout = argout + 1;
  varargout(argout) = {[1200*x(1)^2-400*x(2)+2 -400*x(1); -400*x(1) 200]};
end
return;
