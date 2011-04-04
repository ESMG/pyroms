function [] = viewdata(fname, verbose)
  
    if nargin == 1
        verbose = 1;
    end
  
    if verbose
        fprintf('plotting data from "%s"\n', fname);
        fprintf('  reading %s...', fname);
    end
    data = load(fname);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
    clear data;
  
    if verbose
        fprintf('\n');
    end
  
    if verbose
        fprintf('  plotting...');
    end

    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);
    zmin = min(z);
    zmax = max(z);
    n = length(z);

    map = colormap;
    axis([xmin xmax ymin ymax]);
    axis square;
    set(gca, 'box', 'on');
    hold on;
    for i = 1 : n
        plot(x(i), y(i), 's-', 'color', zcolor(z(i), zmin, zmax, map), 'markersize', 2);
    end
  
    if verbose
        fprintf('\n');
    end

    return
  
function c = zcolor(z, zmin, zmax, map)

    ind = floor((z - zmin) / (zmax - zmin) * 64 + 1);
    ind = min(ind, 64);
    c = map(ind, :);

    return
