BEGIN {
    r = 1;
    T = 10000000;
    N = 100;
}

{
    for (i = 0; i < N; ++i) {
	x = random();
	y = random();
	xx = x * 9.0;
	yy = y * 9.0;
	z = 0.75 * exp(- (xx-2) * (xx-2) / 4 - (yy-2) * (yy-2) / 4) + 0.75 * exp(- (xx-2) * (xx-2) / 49 - (yy-2) / 10) + 0.5 * exp(- (xx-7) * (xx-7) / 4 - (yy-3) * (yy-3) / 4) - 0.2 * exp(- (xx-4) * (xx-4) - (yy-7)*(yy-7));
        printf("%.10g %.10g %.10g\n", x, y, z);
    }
}

END {
}

# One could use in-built generator rand(), but its output may depend on the
# awk implementation used...
#
function random()
{
    r = (r * 40353607) % T;
    return r / T;
}
