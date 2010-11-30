from lpsolve55 import *
from lp_maker import *
f = [110*1.3, 30*2.0, 125*1.56, 75*1.8, 95*.95, 100*2.25, 50*1.35]
A = [[120, 210, 150.75, 115, 186, 140, 85], [110, 30, 125, 75, 95, 100, 50], [1, 1, 1, 1, 1, 1, 1], [1, -1, 0, 0, 0, 0, 0], [0, 0, 1, 0, -2, 0, 0], [0, 0, 0, -1, 0, -1, 1]]
b = [55000, 40000, 400, 0, 0, 0]
lp = lp_maker(f, A, b, [-1, -1, -1, -1, -1, -1], [10, 10, 10, 10, 20, 20, 20], [100, Infinite, 50, Infinite, Infinite, 250, Infinite], None, 1, 0)
solvestat = lpsolve('solve', lp)
obj = lpsolve('get_objective', lp)
print obj
x = lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp', lp)
