function viewexample()
  fprintf(1, '  Drawing results of reconstruction of two data sets with known error:\n');
  fprintf(1, '  (1) data set formed by a superposition of equal number of random points\n');
  fprintf(1, '      of Franke function (std = 0.2) and a linear function (std = 5.0)...');
  subplot(3, 2, 1);
  view('data-F-100.txt', 'out-F-100.txt', 1);
  title('(1): 100+100 data points');
  pause(0.1);
  subplot(3, 2, 3);
  view('data-F-300.txt', 'out-F-300.txt', 1);
  title('(1): 300+300 data points');
  pause(0.1);
  subplot(3, 2, 5);
  view('data-F-1000.txt', 'out-F-1000.txt', 1);
  title('(1): 1000+1000 data points');
  pause(0.1);
  fprintf(1, 'done\n');
  
  fprintf(1, '  (2) data set formed by a superposition of equal number of random points\n');
  fprintf(1, '      of Franke function (std = 5.0) and a linear function (std = 0.2)...');
  subplot(3, 2, 2);
  view('data-P-100.txt', 'out-P-100.txt', 2);
  title('(2): 100+100 data points');
  pause(0.1);
  subplot(3, 2, 4);
  view('data-P-300.txt', 'out-P-300.txt', 2);
  title('(2): 300+300 data points');
  pause(0.1);
  subplot(3, 2, 6);
  view('data-P-1000.txt', 'out-P-1000.txt', 2);
  title('(2): 1000+1000 data points');
  fprintf(1, 'done\n');
  return

function view(data, output, type)
  N = 256;
  for i = 0:20
    V(i+1) = i * 0.1;
  end
  points = load(output);
  k = 1;
  x = zeros(N, N);
  y = zeros(N, N);
  z = zeros(N, N);
  z1 = zeros(N, N);
  for j = N:-1:1
    for i = 1:N
      x(j, i) = points(k, 1);
      y(j, i) = points(k, 2);
      z(j, i) = points(k, 3);
      xx = x(j, i) * 9.0;
      yy = y(j, i) * 9.0;
      z1(j, i) = 0.75 * exp(- (xx-2.0) * (xx-2.0) / 4.0 - (yy-2.0) * (yy-2.0) / 4.0) + 0.75 * exp(- (xx-2.0) * (xx-2.0) / 49.0 - (yy-2.0) / 10.0) + 0.5 * exp(- (xx-7.0) * (xx-7.0) / 4.0 - (yy-3.0) * (yy-3.0) / 4.0) - 0.2 * exp(- (xx-4.0) * (xx-4.0) - (yy-7.0)*(yy-7.0));
      z2(j, i) = x(j, i) + y(j, i);
      k = k + 1;
    end
  end
  
  if (~exist('type', 'var') | type == 0)
    contour(x, y, z1, V, 'k');
    hold on;
  elseif (type == 1)
      contour(x, y, z1, V, 'k');
      hold on;
      contour(x, y, z2, V, 'k:');
      hold on;
  elseif (type == 2)
      contour(x, y, z1, V, 'k:');
      hold on;
      contour(x, y, z2, V, 'k');
      hold on;
  end
  [c, h] = contour(x, y, z, V);
  clabel(c, h);

  points = load(data);
  x = points(:, 1);
  y = points(:, 2);
  plot(x, y, 'k+', 'markersize', 3);
  axis([0 1 0 1]);
  axis square;
  
  return
