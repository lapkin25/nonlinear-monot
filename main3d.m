F1 = @(x,y,z) 3*x + 1.5*y + z - 10;
F2 = @(x,y,z) 0.5*x + 3*y + z - 15;
F3 = @(x,y,z) -50 + 15 * ((x+10).^2 .* sign(x+10) ./100 + (y+10).^2 .* sign(y+10) ./100) + z;

# Решает уравнение F(x0, y) = 0
function y = find_y (F, x0)
  y = solve1dmonot(@(yy) F(x0, yy), 0, 1e-5, 1e-3, 1e5);
endfunction

# Рисует график уравнений F1(x, y) = 0, F2(x, y) = 0, F3(x, y) = 0 на сетке xgrid
function plot_3_2d_functions (F1, F2, F3, xgrid, coord1, coord2)
  fgrid1 = arrayfun(@(x) find_y(F1, x), xgrid);
  fgrid2 = arrayfun(@(x) find_y(F2, x), xgrid);
  fgrid3 = arrayfun(@(x) find_y(F3, x), xgrid);
  figure
  plot(xgrid, fgrid1, xgrid, fgrid2, xgrid, fgrid3)
  xlabel(coord1)
  ylabel(coord2)
endfunction

function plot3dfun (F1, F2, F3)
  xgrid = -20:1:60;
  ygrid = -20:1:60;
  [X, Y] = meshgrid(xgrid, ygrid);
  Z1 = arrayfun(@(x,y) solve1dmonot(@(z) F1(x, y, z), 0, 1e-5, 1e-3, 1e5), X, Y);
  Z2 = arrayfun(@(x,y) solve1dmonot(@(z) F2(x, y, z), 0, 1e-5, 1e-3, 1e5), X, Y);
  Z3 = arrayfun(@(x,y) solve1dmonot(@(z) F3(x, y, z), 0, 1e-5, 1e-3, 1e5), X, Y);
  mesh(X, Y, Z1, 'FaceAlpha', 0.3, 'EdgeColor', 'red')
  hold on
  mesh(X, Y, Z2, 'FaceAlpha', 0.3, 'EdgeColor', 'blue')
  hold on
  mesh(X, Y, Z3, 'FaceAlpha', 0.3, 'EdgeColor', 'green')
  xlabel('x')
  ylabel('y')
  zlabel('z')
endfunction

#plot3dfun(F1, F2, F3);

function y = f (z, F1, F2, F3)
  y = zeros (3, 1);
  y(1) = F1(z(1), z(2), z(3));
  y(2) = F2(z(1), z(2), z(3));
  y(3) = F3(z(1), z(2), z(3));
endfunction

[xx, fval, info] = fsolve (@(z) f(z, F1, F2, F3), [0; 0; 0])
[xx, fval, info] = fsolve (@(z) f(z, F1, F2, F3), [-10; -10; 60])

x0 = 0;
y0 = 0;
z0 = 0;
xgrid = linspace(-60, 60, 50);
ygrid = linspace(-60, 60, 50);
zgrid = linspace(-60, 120, 60);

# manual verification of the 3D algorithm
# choose the following intersection points:
# xy: (1, 2) -> 1  [2]
# yz: (2, 3) -> 2  [2]
# zx: (1, 3) -> 2  [4]
#plot_3_2d_functions(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), @(x,y) F3(x, y, z0), xgrid, "x", "y");
#solve2dmonot(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), x0, y0, 1e-5, 1e-3, 1e5)
x0 = 0.9091;
#plot_3_2d_functions(@(y,z) F1(x0, y, z), @(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), ygrid, "y", "z");
#solve2dmonot(@(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), y0, z0, 1e-5, 1e-3, 1e5)
y0 = -4.1660;
#plot_3_2d_functions(@(z,x) F1(x, y0, z), @(z,x) F2(x, y0, z), @(z,x) F3(x, y0, z), zgrid, "z", "x");
#solve2dmonot(@(z,x) F1(x, y0, z), @(z,x) F3(x, y0, z), z0, x0, 1e-5, 1e-3, 1e5)
z0 = 44.8626;

#plot_3_2d_functions(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), @(x,y) F3(x, y, z0), xgrid, "x", "y");
#solve2dmonot(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), x0, y0, 1e-5, 1e-3, 1e5)
x0 = -7.2477;
#plot_3_2d_functions(@(y,z) F1(x0, y, z), @(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), ygrid, "y", "z");
#solve2dmonot(@(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), y0, z0, 1e-5, 1e-3, 1e5)
y0 = -10.080;
plot_3_2d_functions(@(z,x) F1(x, y0, z), @(z,x) F2(x, y0, z), @(z,x) F3(x, y0, z), zgrid, "z", "x");
solve2dmonot(@(z,x) F1(x, y0, z), @(z,x) F3(x, y0, z), z0, x0, 1e-5, 1e-3, 1e5)
# the root z = 49 is missed!
z0 = 49.4712;

#plot_3_2d_functions(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), @(x,y) F3(x, y, z0), xgrid, "x", "y");
#solve2dmonot(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), x0, y0, 1e-5, 1e-3, 1e5)
x0 = -9.4952;
#plot_3_2d_functions(@(y,z) F1(x0, y, z), @(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), ygrid, "y", "z");
#solve2dmonot(@(y,z) F2(x0, y, z), @(y,z) F3(x0, y, z), y0, z0, 1e-5, 1e-3, 1e5)
y0 = -10.072;
#plot_3_2d_functions(@(z,x) F1(x, y0, z), @(z,x) F2(x, y0, z), @(z,x) F3(x, y0, z), zgrid, "z", "x");
#solve2dmonot(@(z,x) F1(x, y0, z), @(z,x) F3(x, y0, z), z0, x0, 1e-5, 1e-3, 1e5)
z0 = 49.4712;

#plot_3_2d_functions(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), @(x,y) F3(x, y, z0), xgrid, "x", "y");
#solve2dmonot(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), x0, y0, 1e-5, 1e-3, 1e5)

