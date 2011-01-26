from lpsolve55 import *

lp = lpsolve('make_lp', 0, 4)
#lpsolve('set_verbose', lp, IMPORTANT)
lpsolve('set_obj_fn', lp, [1, 3, 6.24, 0.1])
lpsolve('add_constraint', lp, [0, 78.26, 0, 2.9], GE, 92.3)
lpsolve('add_constraint', lp, [0.24, 0, 11.31, 0], LE, 14.8)
lpsolve('add_constraint', lp, [12.68, 0, 0.08, 0.9], GE, 4)
lpsolve('set_lowbo', lp, [28.6, 0, 0, 18])
lpsolve('set_upbo', lp, [Infinite, Infinite, Infinite, 48.98])
lpsolve('set_col_name', lp, ['COLONE', 'COLTWO', 'COLTHREE', 'COLFOUR'])
lpsolve('set_row_name', lp, ['THISROW', 'THATROW', 'LASTROW'])
lpsolve('write_lp', lp, 'a.lp')
print lpsolve('get_mat', lp)[0]
lpsolve('solve', lp)
print lpsolve('get_objective', lp)
print lpsolve('get_variables', lp)[0]
print lpsolve('get_constraints', lp)[0]
lpsolve('delete_lp', lp)
