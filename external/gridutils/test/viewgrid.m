function vg(fname)
  if (~exist(fname, 'file'))
    fprintf(1, '"%s" not found\n', fname);
    return
  end
  
  [nx, ny, x, y] = fload(fname);

  xmin = min(x);
  xmax = max(x);
  ymin = min(y);
  ymax = max(y);
  xdiff = xmax - xmin;
  ydiff = ymax - ymin;
  diff = max([xdiff ydiff]) / 20.0;
  xmin = xmin - diff;
  xmax = xmax + diff;
  ymin = ymin - diff;
  ymax = ymax + diff;
  
  figure;
  hold on;
  axis([xmin xmax ymin ymax]);
  axis image;
  
  fill([xmin xmin xmax xmax], [ymin ymax ymax ymin], [0.7 0.7 0.7]);
  
  ix = 1;
  for (j = 1:ny-1)
    for (i = 1:nx-1)
      i1 = ix;
      i2 = ix + 1;
      i3 = ix + nx + 1;
      i4 = i3 - 1;
      xx = [x(i1) x(i2) x(i3) x(i4)];
      yy = [y(i1) y(i2) y(i3) y(i4)];
      sum = x(i1) + x(i2) + x(i3) + x(i4);
      if ~isnan(sum)
	fill(xx, yy, [1 1 1]);
      end
      ix = ix + 1;
    end
    ix = ix + 1;
  end

  %  plot(x, y, '.');

  return

function [nx, ny, x,y] = fload(fname);
  f = fopen(fname);
  if (f < 0)
    return;
  end
  
  nx = 0;
  ny = 0;
  x = [];
  y = [];
  
  str = fgetl(f);
  if (str >= 0)
    seps = char([32,92]);
    [tmp,str] = strtok(str, seps);
    if (strcmpi(tmp, '##'))
      [nxstr,str] = strtok(str, seps);
      [tmp,str] = strtok(str, seps);
      if (strcmpi(tmp, 'x'))
	nystr = strtok(str, seps);
	nx = sscanf(nxstr, '%d');
	ny = sscanf(nystr, '%d');
      end
    end
  end
  
  fprintf(1, '  reading "%s": %d x %d ... ', fname, nx, ny);

  if (nx > 0 & ny > 0)
    % allocate memory -- is there a better way?
    n = nx * ny;
    x = linspace(0, 0, n);
    y = linspace(0, 0, n);

    str = fgetl(f);
    i = 1;
    while (str >= 0)
      if (str(1) ~= '#')
	xy = sscanf(str, '%f');
	if (length(xy) >= 2)
	  x(i) = xy(1);
	  y(i) = xy(2);
	else
	  x(i) = NaN;
	  y(i) = NaN;
	end
      else
	x(i) = NaN;
	y(i) = NaN;
      end
      i = i + 1;
      if (i > n)
	break;
      end
      str = fgetl(f);
    end
  end
  
  fclose(f);
  fprintf(1, 'ok\n');
  return
