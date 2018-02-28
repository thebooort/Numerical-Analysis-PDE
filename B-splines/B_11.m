function [x] = B_11 (t)
  if (0<=t && t<=1)
      x=(1.-t)^3;
  else
      x=0;
  endif
endfunction