function [x] = B_0 (t)
  if (0<=t && t<=1)
      x=-2*(1.-t)^3+(0.25)*(2-t)^3;
  elseif (1<t && t<=2)
      x=(0.25)*(2-t)^3;
  else
      x=0;
  endif
endfunction