from Numeric import *
from LinearAlgebra import *
x = matrixmultiply(inverse(array([[1, 1], [110, 30]])), array([75, 4000]))
print x
P = matrixmultiply(array([143, 60]), x)
print P