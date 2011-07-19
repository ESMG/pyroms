figure

subplot(2, 2, 1);
viewinterp2('data.txt', 'lin.txt');
title('Linear interpolation (contours)');
pause(0.1);

subplot(2, 2, 2);
viewinterp('data.txt', 'lin.txt');
title('Linear interpolation (colour coded image)');
pause(0.1);

subplot(2, 2, 3);
viewinterp2('data.txt', 'nn.txt');
title(sprintf('NN interpolation (contours)'));
pause(0.1);

subplot(2, 2, 4);
viewinterp('data.txt', 'nn.txt');
title(sprintf('NN interpolation (colour coded image)'));
pause(0.1);

suptitle('Interpolation from elevation contours using  nnbathy');
