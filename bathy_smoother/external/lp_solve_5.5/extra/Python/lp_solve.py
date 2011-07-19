from lpsolve55 import *

def lp_solve(f = None, a = None, b = None, e = None, vlb = None, vub = None, xint = None, scalemode = None, keep = None):
  """LP_SOLVE  Solves mixed integer linear programming problems.

  SYNOPSIS: [obj,x,duals,stat] = lp_solve(f,a,b,e,vlb,vub,xint,scalemode,keep)

     solves the MILP problem

             max v = f'*x
               a*x <> b
                 vlb <= x <= vub
                 x(int) are integer

  ARGUMENTS: The first four arguments are required:

           f: n vector of coefficients for a linear objective function.
           a: m by n matrix representing linear constraints.
           b: m vector of right sides for the inequality constraints.
           e: m vector that determines the sense of the inequalities:
                     e(i) = -1  ==> Less Than
                     e(i) =  0  ==> Equals
                     e(i) =  1  ==> Greater Than
         vlb: n vector of lower bounds. If empty or omitted,
              then the lower bounds are set to zero.
         vub: n vector of upper bounds. May be omitted or empty.
        xint: vector of integer variables. May be omitted or empty.
   scalemode: scale flag. Off when 0 or omitted.
        keep: Flag for keeping the lp problem after it's been solved.
              If omitted, the lp will be deleted when solved.

  OUTPUT: A nonempty output is returned if a solution is found:

         obj: Optimal value of the objective function.
           x: Optimal value of the decision variables.
       duals: solution of the dual problem."""

  if f == None:
          help(lp_solve)
          return

  m = len(a)
  n = len(a[0])
  lp = lpsolve('make_lp', m, n)
  lpsolve('set_verbose', lp, IMPORTANT)
  lpsolve('set_mat', lp, a)
  lpsolve('set_rh_vec', lp, b)
  lpsolve('set_obj_fn', lp, f)
  lpsolve('set_maxim', lp) # default is solving minimum lp.

  for i in range(m):
    if e[i] < 0:
          con_type = LE
    elif e[i] == 0:
          con_type = EQ
    else:
          con_type = GE
    lpsolve('set_constr_type', lp, i + 1, con_type)

  if vlb != None:
    for i in range(n):
      lpsolve('set_lowbo', lp, i + 1, vlb[i])

  if vub != None:
    for i in range(n):
      lpsolve('set_upbo', lp, i + 1, vub[i])

  if xint != None:
    for i in range(len(xint)):
      lpsolve('set_int', lp, xint[i], 1)

  if scalemode != None:
    if scalemode != 0:
      lpsolve('set_scaling', lp, scalemode)

  result = lpsolve('solve', lp)
  if result == 0 or result == 1 or result == 11 or result == 12:
    [obj, x, duals, ret] = lpsolve('get_solution', lp)
    stat = result
  else:
    obj = []
    x = []
    duals = []
    stat = result

  if keep != None and keep != 0:
    lpsolve('delete_lp', lp)

  return [obj, x, duals]
