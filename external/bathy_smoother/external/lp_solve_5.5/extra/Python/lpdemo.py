# Script to demonstrate use of the lp_solve toolkit

from lpsolve55 import *

lp=lpsolve('make_lp',0,4)
lpsolve('add_constraint',lp,[3, 2, 2, 1],3,4)
lpsolve('add_constraint',lp,[0, 4, 3, 1],2,3)
lpsolve('set_obj_fn',lp,[2, 3, -2, 3])
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x

# Change a single element, and maximize

lpsolve('set_mat',lp,2,1,0.5)
lpsolve('set_maxim',lp)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x

# Change RHS

lpsolve('set_rh',lp,1,7.45)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x

# Set var 4 to an integer

lpsolve('set_int',lp,4,1)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x

# Put in lower and upper bounds

lpsolve('set_lowbo',lp,2,2)
lpsolve('set_upbo',lp,4,5.3)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x

# Delete a constraint

lpsolve('del_constraint',lp,1)
lpsolve('add_constraint',lp,[1, 2, 1, 4],3,8)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp',lp)

#%%%%%%%%%%%%

# More examples

# ex1.lp from the lp_solve distribution

lp=lpsolve('make_lp',2,2)
lpsolve('set_mat',lp,[[2, 1],[-4, 4]])
lpsolve('set_obj_fn',lp,[-1, 2])
lpsolve('set_int',lp,[1,1])
lpsolve('set_rh_vec',lp,[5, 5])
lpsolve('set_maxim',lp)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp',lp)

# Example 2

f = [50, 100]
A = [[10, 5],[4, 10], [1, 1.5]]
b = [2500, 2000, 450]
e = [-1, -1, -1]

m = len(A)
n = len(A[0])
lp=lpsolve('make_lp',m,n)
lpsolve('set_obj_fn',lp,f)
lpsolve('set_mat',lp,A)
lpsolve('set_rh_vec',lp,b)
lpsolve('set_maxim',lp)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp',lp)

# Example 3

f = [-40, -36]
vub = [8, 10]
A = [[5, 3]]
b = [45]
e = 1

m = len(A)
n = len(A[0])
lp=lpsolve('make_lp',m,n)
lpsolve('set_obj_fn',lp,f)
lpsolve('set_mat',lp,A)
lpsolve('set_rh_vec',lp,b)
lpsolve('set_constr_type',lp,1,2)
lpsolve('set_upbo',lp,1,8)
lpsolve('set_upbo',lp,2,10)
lpsolve('set_maxim',lp)
result=lpsolve('solve',lp)
obj=lpsolve('get_objective', lp)
print obj
x=lpsolve('get_variables', lp)[0]
print x
lpsolve('delete_lp',lp)
