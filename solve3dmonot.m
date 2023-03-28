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
  # Индекс последовательности приближений: (i_xy, j_xy, k_xy, i_yz, j_yz, k_yz, i_zx, j_zx, k_zx)
  # Вычисление разбивается на множество нитей, каждая соответствует своему индексу
  # Каждая нить обновляет текущее приближение, соответствующее индексу, и может создавать новые нити
  # Для каждой нити в каждом сечении берем k-ю точку пересечения линий {F_i = 0} и {F_j = 0}
  # Обработка отдельного индекса реализована в функции process_index

  global indices  # массив векторов из 9 элементов - индексы
  global ind_cnt  # размер массива индексов
  global xyz  # массив векторов из 3 элементов - текущее приближение для каждого индекса (для каждой нити)
  global resid

  ind_cnt = 0;
  indices = [];
  xyz = [];
  resid = [];

  # Начальные приближения: (i_xy, j_xy, 1, i_yz, j_yz, 1, i_zx, j_zx, 1),
  #   где (i, j) = (1, 2) или (1, 3) или (2, 3)
  for i_xy = 1:3
    for j_xy = i_xy+1:3
      for i_yz = 1:3
        for j_yz = i_yz+1:3
          for i_zx = 1:3
            for j_zx = i_zx+1:3
              ind_cnt += 1;
              indices(ind_cnt, :) = [i_xy, j_xy, 1, i_yz, j_yz, 1, i_zx, j_zx, 1];
              xyz(ind_cnt, :) = [x0, y0, z0];
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor

  for step = 1:num_steps
    current_ind_cnt = ind_cnt;
    processed_ind = zeros(current_ind_cnt);  # обработан ли индекс
    for ind = 1:current_ind_cnt
      process_index(ind, F1, F2, F3, num_steps, x_eps, F_eps, R);
      processed_ind(ind) = true;
    endfor
  endfor

  sol = [];
  for i = 1:ind_cnt
    if (i <= length(resid) && resid(i) < F_eps)
      flag = false;
      for j = 1:i-1
        flag = flag || (max(xyz(i, :) - xyz(j, :)) < x_eps) && resid(j) < F_eps;
      endfor
      # flag == false, если такой же точки не было
      if (!flag)
        sol = [sol; xyz(i, :)];
      endif
    endif
  endfor

endfunction

# Обрабатывает индекс с номером ind
function process_index (ind, F1, F2, F3, num_steps, x_eps, F_eps, R)
  global indices
  global ind_cnt
  global xyz
  global resid

  [i_xy, j_xy, k_xy, i_yz, j_yz, k_yz, i_zx, j_zx, k_zx] = num2cell(indices(ind, :)){:};
  [xx, yy, zz] = num2cell(xyz(ind, :)){:};
  F = {F1, F2, F3};

  # сечение xy
  # пересекаем линии {F_i = 0}, {F_j = 0}; i = i_xy, j = j_xy
  intersection_points = solve2dmonot(@(x,y) F{i_xy}(x, y, zz), @(x,y) F{j_xy}(x, y, zz), xx, yy, x_eps, F_eps, R);
  # обновляем координату x - берем k-ю точку пересечения
  if (k_xy <= rows(intersection_points))
    xx = xyz(ind, 1) = intersection_points(k_xy, 1);
  elseif (rows(intersection_points) > 0)
    xx = xyz(ind, 1) = intersection_points(rows(intersection_points), 1);
  endif

  # пытаемся добавить новые индексы
  for k = 1:rows(intersection_points)
    try_adding_index([i_xy, j_xy, k, i_yz, j_yz, k_yz, i_zx, j_zx, k_zx], [intersection_points(k, 1), yy, zz]);
  endfor

  # сечение yz
  # пересекаем линии {F_i = 0}, {F_j = 0}; i = i_yz, j = j_yz
  intersection_points = solve2dmonot(@(y,z) F{i_yz}(xx, y, z), @(y,z) F{j_yz}(xx, y, z), yy, zz, x_eps, F_eps, R);
  # обновляем координату y - берем k-ю точку пересечения
  if (k_yz <= rows(intersection_points))
    yy = xyz(ind, 2) = intersection_points(k_yz, 1);
  elseif (rows(intersection_points) > 0)
    yy = xyz(ind, 2) = intersection_points(rows(intersection_points), 1);
  endif
  # пытаемся добавить новые индексы
  for k = 1:rows(intersection_points)
    try_adding_index([i_xy, j_xy, k_xy, i_yz, j_yz, k, i_zx, j_zx, k_zx], [xx, intersection_points(k, 1), zz]);
  endfor

  # сечение zx
  # пересекаем линии {F_i = 0}, {F_j = 0}; i = i_zx, j = j_zx
  intersection_points = solve2dmonot(@(z,x) F{i_zx}(x, yy, z), @(z,x) F{j_zx}(x, yy, z), zz, xx, x_eps, F_eps, R);
  # обновляем координату z - берем k-ю точку пересечения
  if (k_zx <= rows(intersection_points))
    zz = xyz(ind, 3) = intersection_points(k_zx, 1);
  elseif (rows(intersection_points) > 0)
    zz = xyz(ind, 3) = intersection_points(rows(intersection_points), 1);
  endif
  # пытаемся добавить новые индексы
  for k = 1:rows(intersection_points)
    try_adding_index([i_xy, j_xy, k_xy, i_yz, j_yz, k_yz, i_zx, j_zx, k], [xx, yy, intersection_points(k, 1)]);
  endfor

  #indices(ind, :)
  #[xx, yy, zz]
  resid(ind) = max([abs(F1(xx, yy, zz)), abs(F2(xx, yy, zz)), abs(F3(xx, yy, zz))]);

endfunction

# Добавляет новый индекс index вместе с соответствующим ему приближением point,
#   если такого индекса нет
function try_adding_index (index, point)
  global indices
  global ind_cnt
  global xyz

  is_index_present = false;
  for i = 1:ind_cnt
    flag = true;
    for s = 1:9
      flag = flag && (indices(i, s) == index(s));
    endfor
    # flag == true, если indices(i, :) == index
    is_index_present = is_index_present || flag;
  endfor
  # если такого индекса нет, добавляем
  if (!is_index_present)
    ind_cnt += 1;
    indices(ind_cnt, :) = index;
    xyz(ind_cnt, :) = point;
  endif
endfunction
