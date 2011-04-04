BEGIN {
    r = 1;
    T = 10000000;
    N = 300;
}

{
    for (i = 0; i < N; ++i) {
	x = (int(random() * 121.0) - 10.0) / 100.0;
	y = (int(random() * 121.0) - 10.0) / 100.0;
	z = 5.0 * x - 3.0 * y;
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
