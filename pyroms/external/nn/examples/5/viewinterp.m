function [h] = viewinterp(fin, fout, verbose)

    if nargin == 2
        verbose = 1;
    end
  
    if verbose
        fprintf('plotting data from "%s" and "%s"\n', fin, fout);
        fprintf('  reading "%s"...', fin);
    end

    data = load(fin);
    xin = data(:, 1);
    yin = data(:, 2);

    if verbose
        fprintf('\n  reading "%s"...', fout);
    end
  
    data = load(fout);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
    clear data;
  
    if verbose
        fprintf('\n');
    end

    if verbose
        fprintf('  working out the grid dimensions...')
    end
    n = length(x);
    if x(2) - x(1) ~= 0 & y(2) - y(1) == 0
        xfirst = 1;
        xinc = x(2) > x(1);
        if xinc
            nx = min(find(diff(x) < 0));
        else
            nx = min(find(diff(x) > 0));
        end
        if mod(n, nx) ~= 0
            error(sprintf('\n  Error: could not work out the grid size, n = %d, nx = %d, n / nx = %f\n', n, nx, n / nx));
        end
        ny = n / nx;
        x = x(1 : nx);
        y = y(1 : nx : n);
        z = reshape(z, nx, ny)';
    elseif x(2) - x(1) == 0 & y(2) - y(1) ~= 0
        xfirst = 0;
        yinc = y(2) > y(1);
        if yinc
            ny = min(find(diff(y) < 0));
        else
            ny = min(find(diff(y) > 0));
        end
        if mod(n, ny) ~= 0
            error(sprintf('\n  Error: could not work out the grid size, n = %d, ny = %d, n / ny = %.3f\n', n, ny, n / ny));
        end
        nx = n / ny;
        y = y(1 : ny);
        x = x(1 : ny : n);
        z = reshape(z, ny, nx);
    else
        error('  Error: not a rectangular grid');
    end
    if verbose
        if xfirst
            fprintf('%d x %d, stored by rows\n', nx, ny);
        else
            fprintf('%d x %d, stored by columns\n', nx, ny);
        end
    end
  
    if verbose
        fprintf('  plotting...');
    end

    h = pcolor_ps(x, y, z);
    zrange = [min(min(z)) max(max(z))];
    caxis(zrange);
    set(h, 'LineStyle', 'none');
    axis equal;
    axis tight;
    hold on;
    plot(xin, yin, 'w.', 'markersize', 1);
    
    if verbose
        fprintf('\n');
    end

    return
  
function [h] = pcolor_ps(x, y, A);

    if nargin == 1
        A = x;
        [n, m] = size(A);
        xx = (0.5 : m + 0.5)';
        yy = (0.5 : n + 0.5)';
    elseif nargin == 3
        n = length(y);
        m = length(x);
        A = reshape(A, n, m); % just in case
        xx = getcorners(x);
        yy = getcorners(y);
    else
        error(sprintf('\n  Error: pcolor_ps(): nargin = %d (expected 1 or 3)\n', nargin));
    end
  
    TMP = zeros(n + 1, m + 1);
    TMP(1 : n, 1 : m) = A;
  
    if nargout == 0
        pcolor(xx, yy, TMP);
    else
        h = pcolor(xx, yy, TMP);
    end
    
    return
    
function [c] = getcorners(x)
  
    n = length(x);
    c = zeros(n + 1, 1);
    c(2 : n) = (x(2 : n) + x(1 : n - 1)) / 2;
    c(1) = 2 * x(1) - c(2);
    c(n + 1) = 2 * x(n) - c(n);

    return
