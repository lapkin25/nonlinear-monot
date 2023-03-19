# Принимает на вход монотонно возрастающую функцию F
# x0 - начальное приближение
# x_eps - допустимая точность x
# F_eps - допустимая точность F
# R - максимальное значение x
function sol = solve1dmonot (F, x0, x_eps, F_eps, R)
  sol = [];

  # вычисляем нижнюю границу для корня
  left = x0;
  x_step = x_eps;
  while (F(left) > 0 && abs(left) < 2 * R)
    left -= x_step;
    x_step *= 2;
  endwhile
  if (F(left) > 0)
    #disp("error")
    return
  endif

  # вычисляем верхнюю границу для корня
  right = x0;
  x_step = x_eps;
  while (F(right) < 0 && abs(right) < 2 * R)
    right += x_step;
    x_step *= 2;
  endwhile
  if (F(right) < 0)
    #disp("error")
    return
  endif

  # F(left) <= 0, F(right) >= 0
  while (right - left > x_eps)
    mid = (left + right) / 2;
    if (F(mid) > 0)
      right = mid;
    else
      left = mid;
    endif
  endwhile
  sol = left;
endfunction

