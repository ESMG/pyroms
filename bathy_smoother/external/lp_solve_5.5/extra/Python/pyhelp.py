import os
import sys

setupfile = file("setup.py",'r')

input = sys.argv[1]

lines = setupfile.readlines()

output = file("setup.py", 'w')

for line in lines:
   if line.startswith("                library_dirs"):
      newline = "                library_dirs=['"+str(input)+"'],\n"
      output.write(newline)
   else:
      output.write(line)

setupfile.close()
output.close()
        
