figure

subplot(3, 2, 1);
fname = 'data.txt';
fprintf('plotting data points from "%s"\n', fname);
fprintf('  reading %s...', fname);
data = load(fname);
xrange = [min(data(:, 1)) max(data(:, 1))];
yrange = [min(data(:, 2)) max(data(:, 2))];
zrange = [min(data(:, 3)) max(data(:, 3))];
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
viewinterp('data.txt', 'lin.txt');
title('Linear interpolation');
pause(0.1);

subplot(3, 2, 4);
viewinterp('data.txt', 'nn-inf.txt');
caxis(zrange);
title(sprintf('Natural Neighbours interpolation\n(extrapolation allowed)'));
pause(0.1);

subplot(3, 2, 5);
viewinterp('data.txt', 'nn-0.txt');
caxis(zrange);
title(sprintf('Natural Neighbours interpolation\n(interpolation only)'));
pause(0.1);

subplot(3, 2, 6);
viewinterp('data.txt', 'nn-ns.txt');
caxis(zrange);
title(sprintf('Non-Sibsonian NN interpolation\n(interpolation only)'));
pause(0.1);

suptitle('Interpolation of bathymetry data from sonar using nnbathy');
