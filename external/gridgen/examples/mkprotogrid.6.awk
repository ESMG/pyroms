BEGIN {
# conformal modulus
    mod = 1
    nx = 51;
# maximal aspect ratio
    a = 2.0;
    PI = 3.14159265357989;
# step by X (uniform)
    dx = 1.0 / (nx - 1);
# initial step by Y
    dy0 = dx * mod;
}

{
    y = 0.0;
    for (j = 1; y <= 1.0; ++j) {
	x = 0.0;
	for (i = 1; i <= nx; ++i) {
	    printf("%.10g %.10g\n", x, y);
	    x += dx;
	}
	s = sin(PI * y);
	dy = dy0 * (1.0 + (a - 1.0) * s * s);
	y += dy;
    }

    printf("## %d x %d\n", nx, j - 1);
}
