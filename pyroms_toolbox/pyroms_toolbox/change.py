import numpy as np

def change(old,relation,flag,value):
    """ 
    Change values in a matrix
    """

    if relation != '==' and relation != '!=' and relation != '>' and relation != '<' and \
         relation != '!=' and relation != '>=' and relation != '<=':
        raise ValueError('Relation {%s} not valid' % relation)

    if np.isnan(flag):
        if relation == '==':
            replace = np.where(np.isnan(old) == True)
        elif relation == '!=':
            replace = np.where(np.isnan(old) == False)
        else:
            raise ValueError('Relation should be == or ~= to compare to NaN')

    else:
        if relation == '==':
            replace = np.where(old == flag)
        elif relation == '>':
            replace = np.where(old > flag)
        elif relation == '<':
            replace = np.where(old < flag)
        elif relation == '!=':
            replace = np.where(old != flag)
        elif relation == '>=':
            replace = np.where(old >= flag)
        elif relation == '<=':
            replace = np.where(old <= flag)

    new = old.copy()
    new[replace] = value

    return new
