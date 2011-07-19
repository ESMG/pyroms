from lp_maker import *
f = [143, 60, 195]
A = [[120, 210, 150.75], [110, 30, 125], [1, 1, 1]]
b = [15000, 4000, 75]
lp = lp_maker(f, A, b, [-1, -1, -1], None, None, None, 1, 0)
solvestat = lpsolve('solve', lp)
obj = lpsolve('get_objective', lp)
print obj
x = lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp', lp)
