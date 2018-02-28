function [x] = B_1 (t)
  if (0<=t && t<=1)
      x=(1/6.)*(3-t)^3-(0.75)*(2-t)^3+(1.5)*(1-t)^3;
  elseif (1<t && t<=2)
      x=(1/6.)*(3-t)^3-(0.75)*(2-t)^3;
  elseif (2<t && t<=3)
      x=(1/6.)*(3-t)^3;
  else
      x=0;
  endif
endfunction