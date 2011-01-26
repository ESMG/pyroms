# Test and Demonstrate the use of lp_solve

from lp_solve import *

# Example 1 from the lp_solve distribution
f = [-1, 2]
A = [[2, 1], [-4, 4]]
b = [5, 5]
e = [-1, -1]
xint = [1, 2]

[v,x,duals] = lp_solve(f,A,b,e,None,None,xint)

print v
print x

# Example 2
f = [50, 100]
A = [[10, 5],[4, 10],[1, 1.5]]
b = [2500, 2000, 450]
e = [-1, -1, -1]

[v,x,duals] = lp_solve(f,A,b,e)
print v
print x

# Example 3
f = [-40, -36]
vub = [8, 10]
A = [[5, 3]]
b = [45]
e = [1]

[v,x,duals] = lp_solve(f,A,b,e,None,vub)
print v
print x

# Example 4
f = [10, 6, 4]
A = [[1, 1, 1], [10, 4, 5], [2, 2, 6]]
b = [100, 600, 300]
e = [-1, -1, -1]
xint = [2]

[v,x,duals] = lp_solve(f,A,b,e,None,None,xint)
print v
print x

# Example 5
# Integer programming example, page 218 of Ecker & Kupferschmid
f = [-3, 7, 12]
a = [[-3, 6, 8], [6, -3, 7], [-6, 3, 3]]
b = [12, 8, 5]
e = [-1, -1, -1]
xint = [1, 2, 3]

[v,x,duals] = lp_solve(f,a,b,e,None,None,xint)
print v
print x

# Example 6
# 0-1 programming example, page 228 233 of Ecker & Kupferschmid
f = [-2, -3, -7, -7]
a = [[1, 1, -2, -5], [-1, 2, 1, 4]]
b = [2, -3]
e = [1, 1]
xint = [1, 2, 3, 4]
vub = [1, 1, 1, 1]

[v,x,duals] = lp_solve(f,a,b,e,None,vub,xint)
print v
print x

# Example 7
# 0-1 programming example, page 238 of Ecker & Kupferschmid
f = [-1, -2, -3, -7, -8, -8]
a = [[5, -3, 2, -3, -1, 2], [-1, 0, 2, 1, 3, -3], [1, 2, -1, 0, 5, -1]]
b = [-5, -1, 3]
e = [1, 1, 1]
xint = [1, 2, 3, 4, 5, 6]
vub = [1, 1, 1, 1, 1, 1]

[v,x,duals] = lp_solve(f,a,b,e,None,vub,xint)
print v
print x

# ex2.lp from the lp_solve distribution
f=[8, 15]
a = [[10, 21], [2, 1]]
b = [156, 22]
e = [-1, -1]
[v,x,duals] = lp_solve(f,a,b,e)
print v
print x

# ex3.lp from the lp_solve distribution
f=[3, 13]
a = [[2, 9], [11, -8]]
b = [40, 82]
e = [-1, -1]
[v,x,duals] = lp_solve(f,a,b,e)
print v
print x

# ex6.lp from the lp_solve distribution
f=[592, 381, 273, 55, 48, 37, 23]
a = [[3534, 2356, 1767, 589, 528, 451, 304]]
b = [119567]
e = [-1]
xint = [1, 2, 3, 4, 5, 6, 7]
vub = None
[v,x,duals] = lp_solve(f,a,b,e,None,vub,xint)
print v
print x
