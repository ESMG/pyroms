subplot(2, 2, 1)
vb('bathy-cs.txt','x.txt','y.txt',102,140,60);
title('Approximation by C_3^1 spline');
subplot(2, 2, 2)
vb('bathy-nn.txt','x.txt','y.txt',102,140,60);
title('Natural Neighbours interpolation');
subplot(2, 2, 3)
vb('bathy-ns.txt','x.txt','y.txt',102,140,60);
title('Non-Sibsoninan NN interpolation');
subplot(2, 2, 4)
vb('bathy-l.txt','x.txt','y.txt',102,140,60);
title('Linear interpolation');
