function display(S)
  disp(sprintf('Size: (%d,%d)',S.n(1),S.n(2)));
  [ir,ic] = ind2sub2(S.n,S.i);
  disp('Row  Col     Val')
  for i=1:length(ir)
    disp(sprintf(' %-5d%-5d%.4f', ir(i),ic(i), S.v(i)));    
  end
  
