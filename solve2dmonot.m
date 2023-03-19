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

  # пускай начальное приближение принадлежит линии {F1(x0, y0) = 0}
  x00 = solve1dmonot(@(x) F1(x, y0), x0, x_eps, F_eps, R);
  y00 = y0;
  # вдруг мы попали в точности на корень
  root_twice = is_root(F1, F2, x00, y00, x_eps);

  # пусть x возрастает, y убывает
  x = x00;
  y = y00;
  while (abs(x) < R && abs(y) < R)
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

    # идем дальше, чтобы перейти к поиску следующего корня
    while (is_root(F1, F2, x, y, x_eps))
      x += x_eps;
      y -= x_eps;
    endwhile
  endwhile

  # теперь пусть x убывает, y возрастает
  x = x00;
  y = y00;
  root_twice_flag = !root_twice;
  while (abs(x) < R && abs(y) < R)
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
      if (root_twice_flag)
        sol = [ x, y ; sol ];
      else
        root_twice_flag = true;
      endif
    endif

    # идем дальше, чтобы перейти к поиску следующего корня
    while (is_root(F1, F2, x, y, x_eps))
      x -= x_eps;
      y += x_eps;
    endwhile
  endwhile
endfunction

# проверяет, что x и y находятся вблизи кривых {F1 = 0} и {F2 = 0}
function ret = is_root(F1, F2, x, y, x_eps)
  s1 = sign(F1(x, y));
  s2 = sign(F2(x, y));
  is_root_x = ( (sign(F1(x - x_eps, y)) != s1) || (sign(F1(x + x_eps, y)) != s1) ) ...
    && ( (sign(F2(x - x_eps, y)) != s2) || (sign(F2(x + x_eps, y)) != s2) );
  is_root_y = ( (sign(F1(x, y - x_eps)) != s1) || (sign(F1(x, y + x_eps)) != s1) ) ...
    && ( (sign(F2(x, y - x_eps)) != s2) || (sign(F2(x, y + x_eps)) != s2) );
  ret = is_root_x || is_root_y;
endfunction
