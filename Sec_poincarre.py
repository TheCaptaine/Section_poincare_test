import numpy as np
import sys

data = np.loadtxt(sys.argv[1]).T

"""with open("SP.dat", "w") as file:
	file.write("# t q2 p2\n")
	for i in range(len(data[0])-1):
		try:
			a = float(str(data[1][i])[:3])
		except ValueError:
			#print(1, str(data[1][i])[:3], str(data[1][i+1])[:3])
			try:
				a = float(str(data[1][i])[:2])
			except ValueError:
				#print(2, str(data[1][i])[:2], str(data[1][i+1])[:2])
				a = float(str(data[1][i])[:1])
		try:
			b = float(str(data[1][i+1])[:3])
		except ValueError:
			#print(1, str(data[1][i])[:3], str(data[1][i+1])[:3])
			try:
				b = float(str(data[1][i+1])[:2])
			except ValueError:
				#print(2, str(data[1][i])[:2], str(data[1][i+1])[:2])
				b = float(str(data[1][i+1])[:1])
		if a*b < 0:
			file.write("{} {} {}\n".format(data[0][i], data[2][i]-a*(data[2][i+1]-data[2][i])/(b-a), data[4][i]-a*(data[4][i+1]-data[4][i])/(b-a)))
"""

with open(sys.argv[2], "w") as file:
	file.write("# t q2 p2\n")
	for i in range(len(data[0])-1):
		a, b = data[1][i], data[1][i+1]
		#print(a, b, a*b)
		if a*b < 0:
			file.write("{} {} {}\n".format(data[0][i], data[2][i]-a*(data[2][i+1]-data[2][i])/(b-a), data[4][i]-a*(data[4][i+1]-data[4][i])/(b-a)))

