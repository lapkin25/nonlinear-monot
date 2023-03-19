F1 = @(x,y,z) 3*x + 1.5*y + z - 10;
F2 = @(x,y,z) 0.5*x + 3*y + z - 15;
F3 = @(x,y,z) -50 + 15 * ((x+10).^2 .* sign(x+10) ./100 + (y+10).^2 .* sign(y+10) ./100) + z;

x0 = 0;
y0 = 0;
z0 = 0;

global xyz
global indices
global ind_cnt
xyz = [];
indices = [];
ind_cnt = 0;


solve3dmonot(F1, F2, F3, x0, y0, z0, 2, 1e-3, 1e-3, 1e5)

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

xgrid = linspace(-25, 25, 50);

plot_3_2d_functions(@(x,y) F1(x, y, z0), @(x,y) F2(x, y, z0), @(x,y) F3(x, y, z0), xgrid, "x", "y");

solve2dmonot(@(x,y) F1(x, y, z0), @(x,y) F3(x, y, z0), x0, y0, 1e-3, 1e-3, 1e5)
