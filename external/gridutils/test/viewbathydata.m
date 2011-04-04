function viewbathydata(fname, zmax, txt);
  
  p = load(fname);
  doprint = 0;
  if exist('txt', 'var')
    if txt ~= 0
      doprint = 1;
    end
  end

  if exist('zmax', 'var')
    p = p(find(p(:,3) < zmax),:);
  else
    zmax = max(p(:,3));
  end
  
  zmin = min(0.0, min(p(:,3)));
  if (zmax == zmin)
    zmax = zmin + 1.0;
  end

  hold on;
  
  n = length(p);
  map = colormap;
  if (~doprint)
    for i = 1:n
      pp = p(i,:);
      plot(pp(1), pp(2), 's-', 'color', map(cindex(pp(3), zmin, zmax), :), 'markersize', 2);
    end
  else
    for i = 1:n
      pp = p(i,:);
      plot(pp(1), pp(2), 's-', 'color', map(cindex(pp(3), zmin, zmax), :), 'markersize', 2);
      s = sprintf('%.1f', pp(3));
      text(pp(1), pp(2), s);
    end
  end
  hold off;
  
  title(sprintf('%s', fname));
  h = colorbar;
  set(h, 'YTick', [1 65]);
  set(h, 'YTickLabel', [zmin zmax]);

  return;
  
function index = cindex(z, zmin, zmax)
  
  if z <= zmin
    index = 1;
  elseif z >= zmax
    index = 64;
  else
    index = floor((z - zmin) / (zmax - zmin) * 64 + 1);
  end

  return
