import numpy as np

def np2latex(A,name=None, rounding=3):
	nrows, ncols = A.shape
	print("\\begin{equation}")
	if name != None:
		print("  \mathbf{"+name+"}=",end='')
	print("\\left[\\begin{array}{",end='')
	for c in range(ncols):
		print("c",end='')
	print("}")
	for i in range(nrows):
		for j in range(ncols):
			if j == ncols-1:
			 print("    "+str(round(A[i,j],rounding))+" "+"\\\\")
			elif j == 0:
				print("    "+str(round(A[i,j],rounding))+" "+"&",end='')
			else:
				print("    "+str(round(A[i,j],rounding))+" "+"&",end='')
	print("  \\end{array}\\right]")
	print("\\end{equation}")


