solve2dmonot(@(x,y) x + y + 1, @(x,y) 2*x + y - 1, 0, 0, 1e-5, 1e-3, 1e6)

f1 = @(x,y) abs(x)*x + abs(y)*y - 1;
f2 = @(x,y) x + y - 1;
solve2dmonot(f1, f2, 0, 0, 1e-5, 1e-3, 1000)

f1 = @(x,y) abs(x)*x + abs(y)*y - 1;
f2 = @(x,y) x + y - 0.5;
solve2dmonot(f1, f2, 0, 0, 1e-5, 1e-3, 1000)

