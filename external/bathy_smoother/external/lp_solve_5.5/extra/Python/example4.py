from lp_maker import *
f = [143, 60]
A = [[120, 210], [110, 30], [1, 1]]
b = [15000, 4000, 75]
lp = lp_maker(f, A, b, [-1, -1, -1], None, None, None, 1, 0)
solvestat = lpsolve('solve', lp)
obj = lpsolve('get_objective', lp)
print obj
x = lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp', lp)
