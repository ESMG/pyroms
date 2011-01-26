figure

subplot(3, 2, 1);
fname = 'data.txt';
fprintf('plotting data points from "%s"\n', fname);
fprintf('  reading %s...', fname);
data = load(fname);
xrange = range(data(:, 1));
yrange = range(data(:, 2));
zrange = range(data(:, 3));
fprintf('\n');
fprintf('  plotting...');
axis([xrange yrange]);
axis tight;
axis square;
set(gca, 'box', 'on');
hold on;
plot(data(:, 1), data(:, 2), 'k.', 'markersize', 1);
fprintf('\n');
clear data;
title('Data points');
pause(0.1);

subplot(3, 2, 2);
viewdata('data.txt');
title('Data');
pause(0.1);

subplot(3, 2, 3);
viewinterp('data.txt', 'out_default.txt');
caxis(zrange);
title(sprintf('Approximation with Bivariate Cubic Spline\n(default settings)'));
pause(0.1);

subplot(3, 2, 4);
viewinterp('data.txt', 'out_nppc20.txt');
caxis(zrange);
title(sprintf('Approximation with Bivariate Cubic Spline\n(with -P nppc=20)'));
pause(0.1);

subplot(3, 2, 5);
viewinterp('data.txt', 'out_k70.txt');
caxis(zrange);
title(sprintf('Approximation with Bivariate Cubic Spline\n(with -P k=70)'));
pause(0.1);

subplot(3, 2, 6);
viewinterp('data.txt', 'out_k280.txt');
caxis(zrange);
title(sprintf('Approximation with Bivariate Cubic Spline\n(with -P k=280)'));
pause(0.1);

suptitle('Interpolation of bathymetry data from sonar using csabathy');
