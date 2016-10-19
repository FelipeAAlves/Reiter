function [ir,ic] = ind2sub2(nn,i)
  i = i-1;
  ic = floor(i/nn(1));
  ir = i - ic*nn(1) + 1;
  ic = ic+1;