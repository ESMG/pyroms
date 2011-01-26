function viewbound(file);
  [x, y] = fload(file);

  plot(x,y,'k.-');
  axis equal;
  hold on;

  return;

function [x,y] = fload(fname);
  f = fopen(fname);
  if (f < 0)
    return;
  end
  
  x = [];
  y = [];
  
  str = fgetl(f);
  while (str >= 0)
    if (str(1) ~= '#')
      xy = sscanf(str, '%f');
      if (length(xy) == 2)
	x = [x xy(1)];
	y = [y xy(2)];
      end
    end
    str = fgetl(f);
  end
  
  fclose(f);
  return
 