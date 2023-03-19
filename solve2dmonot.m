# Принимает две функции F1(x, y) и F2(x, y),
#   монотонно возрастающие по x и y,
#   решает систему F1(x, y) = 0, F2(x, y) = 0,
# x0, y0 - начальное приближение
# x_eps - допустимая точность x, y
# F_eps - допустимая точность F
# R - максимальное значение |x| и |y|
# Возвращает матрицу N_roots x 2 из корней системы
#   в порядке возрастания x
function sol = solve2dmonot (F1, F2, x0, y0, x_eps, F_eps, R)
  sol = [];

  # пускай приближение принадлежит линии {F1(x, y) = 0}
  x00 = solve1dmonot(@(xx) F1(xx, y0), x0, x_eps, F_eps, R);
  y00 = y0;
  # вдруг мы попали в точности на корень
  [x00, y00] = skip_root_right(F1, F2, x00, y00, 2 * x_eps, F_eps, R);

  # пусть x возрастает, y убывает
  x = x00;
  y = y00;
  while (abs(x) < R && abs(y) < R)
    # идем дальше, чтобы перейти к поиску следующего корня
    [x, y] = skip_root_right(F1, F2, x, y, 2 * x_eps, F_eps, R);
    # приближение принадлежит линии {F1(x, y) = 0}
    while (abs(x) < R && abs(y) < R && !is_root(F1, F2, x, y, x_eps))
      # F1(x, y) == 0
      if (F2(x, y) > 0)  # y будет убывать
        y = solve1dmonot(@(yy) F2(x, yy), y, x_eps, F_eps, R);
        # F2(x, y) == 0
      else  # x будет возрастать
        x = solve1dmonot(@(xx) F2(xx, y), x, x_eps, F_eps, R);
        # F2(x, y) == 0
      endif

      # F2(x, y) == 0
      if (F1(x, y) > 0)  # y будет убывать
        y = solve1dmonot(@(yy) F1(x, yy), y, x_eps, F_eps, R);
        # F1(x, y) = 0
      else  # x будет возрастать
        x = solve1dmonot(@(xx) F1(xx, y), x, x_eps, F_eps, R);
        # F1(x, y) == 0
      endif
    endwhile

    if (is_root(F1, F2, x, y, x_eps))
      # добавляем найденную пару (x, y) к множеству корней
      sol = [ sol ; x, y ];
    endif
  endwhile

  # теперь пусть x убывает, y возрастает
  x = x00;
  y = y00;
  while (abs(x) < R && abs(y) < R)
    # идем дальше, чтобы перейти к поиску следующего корня
    [x, y] = skip_root_left(F1, F2, x, y, 2 * x_eps, F_eps, R);
    # приближение принадлежит линии {F1(x, y) = 0}
    while (abs(x) < R && abs(y) < R && !is_root(F1, F2, x, y, x_eps))
      # F1(x, y) == 0
      if (F2(x, y) < 0)  # y будет возрастать
        y = solve1dmonot(@(yy) F2(x, yy), y, x_eps, F_eps, R);
        # F2(x, y) == 0
      else  # x будет убывать
        x = solve1dmonot(@(xx) F2(xx, y), x, x_eps, F_eps, R);
        # F2(x, y) == 0
      endif

      # F2(x, y) == 0
      if (F1(x, y) < 0)  # y будет возрастать
        y = solve1dmonot(@(yy) F1(x, yy), y, x_eps, F_eps, R);
        # F1(x, y) = 0
      else  # x будет убывать
        x = solve1dmonot(@(xx) F1(xx, y), x, x_eps, F_eps, R);
        # F1(x, y) == 0
      endif
    endwhile

    if (is_root(F1, F2, x, y, x_eps))
      # добавляем найденную пару (x, y) к множеству корней
      sol = [ x, y ; sol ];
    endif
  endwhile
endfunction

function [x, y] = skip_root_right (F1, F2, x0, y0, x_eps, F_eps, R)
  x = x0;
  y = y0;
  while (is_root(F1, F2, x, y, x_eps))
    while (is_root(F1, F2, x, y, x_eps))
      x += x_eps;
      y -= x_eps;
    endwhile
    # пускай приближение принадлежит линии {F1(x, y) = 0}
    x = solve1dmonot(@(xx) F1(xx, y), x, x_eps, F_eps, R);
  endwhile
  # теперь is_root(F1, F2, x, y, x_eps) == false
endfunction

function [x, y] = skip_root_left (F1, F2, x0, y0, x_eps, F_eps, R)
  x = x0;
  y = y0;
  while (is_root(F1, F2, x, y, x_eps))
    while (is_root(F1, F2, x, y, x_eps))
      x -= x_eps;
      y += x_eps;
    endwhile
    # пускай приближение принадлежит линии {F1(x, y) = 0}
    x = solve1dmonot(@(xx) F1(xx, y), x, x_eps, F_eps, R);
  endwhile
  # теперь is_root(F1, F2, x, y, x_eps) == false
endfunction

# проверяет, что x и y находятся вблизи кривых {F1 = 0} и {F2 = 0}
function ret = is_root (F1, F2, x, y, x_eps)
  s1 = sign(F1(x, y));
  s2 = sign(F2(x, y));
  is_root_F1 = (sign(F1(x - x_eps, y)) != s1) || (sign(F1(x + x_eps, y)) != s1) ...
    || (sign(F1(x, y - x_eps)) != s1) || (sign(F1(x, y + x_eps)) != s1);
  is_root_F2 = (sign(F2(x - x_eps, y)) != s2) || (sign(F2(x + x_eps, y)) != s2) ...
    || (sign(F2(x, y - x_eps)) != s2) || (sign(F2(x, y + x_eps)) != s2);
  ret = is_root_F1 && is_root_F2;
endfunction
