x = randn(5,10);
y = ssparse(randn(10,3));

z = x*y;
z2 = mtimes_old(x,y);
distance(full(z),full(z2));

