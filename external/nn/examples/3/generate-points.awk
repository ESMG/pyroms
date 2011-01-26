BEGIN {
    N = 101;
    inc = 1.0 / (N - 1);
}

{
    y = 0.0;
    for (j = 0; j < N; ++j) {
	x = 0.0;
	for (i = 0; i < N; ++i) {
	    printf("%.2f %.2f %.2f\n", x, y, 5.0 * x - 3.0 * y);
	    x += inc;
	}
	y += inc;
    }
}

END {
}
