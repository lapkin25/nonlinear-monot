# Принимает три функции F1(x, y, z) и F2(x, y, z), F3(x, y, z)
#   монотонно возрастающие по x, y, z,
#   решает систему F1(x, y, z) = 0, F2(x, y, z) = 0, F3(x, y, z) = 0
# x0, y0, z0 - начальное приближение
# num_steps - число шагов алгоритма
# x_eps - допустимая точность x, y, z
# F_eps - допустимая точность F
# R - максимальное значение |x|, |y|, |z|
# Возвращает матрицу N_roots x 3 из корней системы
function sol = solve3dmonot (F1, F2, F3, x0, y0, z0, num_steps, x_eps, F_eps, R)
  global indices
  global ind_cnt
  global xyz

  sol = [];

  # множество индексов вида [i, j, k] --
  #   "по x выбрать i-ю точку пересечения пар линий нулевого уровня, по y - j-ю, по z - k-ю"
  indices = [];
  ind_cnt = 0;

  # текущие приближения для каждого такого индекса
  xyz = [];

  # добавляем начальное приближение к xyz (индексы 1, 1, 1)
  ind_cnt += 1;
  indices(ind_cnt, :) = [1, 1, 1];
  xyz(ind_cnt, :) = [x0, y0, z0];

  # debug
  ind_cnt += 1;
  indices(ind_cnt, :) = [1, 3, 3];
  xyz(ind_cnt, :) = [x0, y0, z0];

  for step = 1:num_steps
    current_ind_cnt = ind_cnt;
    processed_ind = zeros(current_ind_cnt);  # обработан ли индекс
    resid_ind = zeros(current_ind_cnt);  # невязки текущих приближений
    for ind = 1:current_ind_cnt
      [xx, yy, zz] = num2cell(xyz(ind, :)){:};
      resid_ind(ind) = max([abs(F1(xx, yy, zz)), abs(F2(xx, yy, zz)), abs(F3(xx, yy, zz))]);
    endfor
    while true
      # сортировка выбором: выбираем на каждом шаге индекс с минимальной невязкой
      best_ind_set = false;
      for ind = 1:current_ind_cnt
        if (!processed_ind(ind) && (!best_ind_set || resid_ind(best_ind) > resid_ind(ind)))
          best_ind = ind;
          best_ind_set = true;
        endif
      endfor
      if (!best_ind_set)  # необработанные индексы кончились
        break
      endif
      ind = best_ind;
      processed_ind(ind) = true;
      process_index(ind, F1, F2, F3, num_steps, x_eps, F_eps, R);

    endwhile
  endfor

endfunction

# Обрабатывает индекс с номером ind
function process_index (ind, F1, F2, F3, num_steps, x_eps, F_eps, R)
  global indices
  global ind_cnt
  global xyz

  [xx, yy, zz] = num2cell(xyz(ind, :)){:};
  [xi, yi, zi] = num2cell(indices(ind, :)){:};

  for step = 1:num_steps
    # сечение xy
    # попарно пересекаем линии {F1 = 0}, {F2 = 0}, {F3 = 0}
    pairwise_intersection_points = [solve2dmonot(@(x,y) F1(x, y, zz), @(x,y) F2(x, y, zz), xx, yy, x_eps, F_eps, R) ; ...
      solve2dmonot(@(x,y) F1(x, y, zz), @(x,y) F3(x, y, zz), xx, yy, x_eps, F_eps, R) ; ...
      solve2dmonot(@(x,y) F2(x, y, zz), @(x,y) F3(x, y, zz), xx, yy, x_eps, F_eps, R) ];
    pairwise_intersection_points = ...
      [sort(pairwise_intersection_points(:, 1)), sort(pairwise_intersection_points(:, 2), "descend")];
    # обновляем координату x у приближения с текущим индексом
    if (xi <= length(pairwise_intersection_points))
      xx = xyz(ind, 1) = pairwise_intersection_points(xi);
    endif
    # пытаемся добавить новые индексы
    for i = 1:rows(pairwise_intersection_points)
      try_adding_index([i, yi, zi], [pairwise_intersection_points(i, 1), yy, zz]);
    endfor

    # сечение yz
    # попарно пересекаем линии {F1 = 0}, {F2 = 0}, {F3 = 0}
    pairwise_intersection_points = [solve2dmonot(@(y,z) F1(xx, y, z), @(y,z) F2(xx, y, z), yy, zz, x_eps, F_eps, R) ; ...
      solve2dmonot(@(y,z) F1(xx, y, z), @(y,z) F3(xx, y, z), yy, zz, x_eps, F_eps, R) ; ...
      solve2dmonot(@(y,z) F2(xx, y, z), @(y,z) F3(xx, y, z), yy, zz, x_eps, F_eps, R) ];
    pairwise_intersection_points = ...
      [sort(pairwise_intersection_points(:, 1)), sort(pairwise_intersection_points(:, 2), "descend")];
    # обновляем координату y у приближения с текущим индексом
    if (yi <= length(pairwise_intersection_points))
      yy = xyz(ind, 2) = pairwise_intersection_points(yi);
    endif
    # пытаемся добавить новые индексы
    for i = 1:rows(pairwise_intersection_points)
      try_adding_index([xi, i, zi], [xx, pairwise_intersection_points(i, 1), zz]);
    endfor

    # сечение zx
    # попарно пересекаем линии {F1 = 0}, {F2 = 0}, {F3 = 0}
    pairwise_intersection_points = [solve2dmonot(@(z,x) F1(x, yy, z), @(z,x) F2(x, yy, z), zz, xx, x_eps, F_eps, R) ; ...
      solve2dmonot(@(z,x) F1(x, yy, z), @(z,x) F3(x, yy, z), zz, xx, x_eps, F_eps, R) ; ...
      solve2dmonot(@(z,x) F2(x, yy, z), @(z,x) F3(x, yy, z), zz, xx, x_eps, F_eps, R) ];
    pairwise_intersection_points = ...
      [sort(pairwise_intersection_points(:, 1)), sort(pairwise_intersection_points(:, 2), "descend")];
    # обновляем координату z у приближения с текущим индексом
    if (zi <= length(pairwise_intersection_points))
      zz = xyz(ind, 3) = pairwise_intersection_points(zi);
    endif
    # пытаемся добавить новые индексы
    for i = 1:rows(pairwise_intersection_points)
      try_adding_index([xi, yi, i], [xx, yy, pairwise_intersection_points(i, 1)]);
    endfor

    #step
    #xyz

    #[xx, yy, zz]
    #resid = max([abs(F1(xx, yy, zz)), abs(F2(xx, yy, zz)), abs(F3(xx, yy, zz))])
    #if (resid > resid_ind(ind))
    #  disabled_ind(ind) = true;
    #endif
  endfor

  indices(ind, :)
  [xx, yy, zz]
  resid = max([abs(F1(xx, yy, zz)), abs(F2(xx, yy, zz)), abs(F3(xx, yy, zz))])

endfunction

# Добавляет новый индекс index вместе с соответствующим ему приближением point,
#   если такого индекса нет
function try_adding_index (index, point)
  global indices
  global ind_cnt
  global xyz

  is_index_present = false;
  for i = 1:ind_cnt
    if (indices(i, 1) == index(1) && indices(i, 2) == index(2) && indices(i, 3) == index(3))
      is_index_present = true;
    endif
  endfor
  # если такого индекса нет, добавляем
  if (!is_index_present)
    ind_cnt += 1;
    indices(ind_cnt, :) = index;
    xyz(ind_cnt, :) = point;
  endif
endfunction
