function [x2,p2] = ptd_percrange(x,p,range)

sump = cumsum([0;p]);
sump = sump./sump(end);

ilo = lookup(sump,range(1),0);
droplo = (range(1)-sump(ilo))/p(ilo);
incllo = 1 - droplo;

n = length(p);
if(range(2)==1)
  ihi = n;
  inclhi = 1;
else
  ihi = lookup(sump,range(2),0);
  inclhi = (range(2)-sump(ihi))/p(ihi);
  if(inclhi==0)
    ihi = ihi-1;
    inclhi = 1;
  end
end

x2 = x(ilo:ihi);
p2 = p(ilo:ihi);
if(ilo==ihi)
  p2 = p2*(incllo+inclhi-1);
else
  p2(1) = p2(1)*incllo;
  p2(end) = p2(end)*inclhi;
end

assert(abs(diff(range)-sum(p2))<1e-12);




