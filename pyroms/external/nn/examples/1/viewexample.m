function viewexample()

    subplot(3, 1, 1);
    view('data-100.txt', 'out-100.txt');
    title('100 data points');
    pause(0.1);

    subplot(3, 1, 2);
    view('data-300.txt', 'out-300.txt');
    title('300 data points');
    pause(0.1);

    subplot(3, 1, 3);
    view('data-1000.txt', 'out-1000.txt');
    title('1000 data points');

    suptitle('Interpolation of Franke test function using nnbathy');

    return

function view(data, output)

    N = 256;
    for i = 0:15
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
            k = k + 1;
        end
    end
  
    contour(x, y, z1, V, 'k');
    hold on;
    [c, h] = contour(x, y, z, V);
    clabel(c, h);

    points = load(data);
    x = points(:, 1);
    y = points(:, 2);
    plot(x, y, 'k+', 'markersize', 3);
    axis([0 1 0 1]);
    axis square;
  
    return
