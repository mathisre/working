import os
import sys
import numpy as np

try:
	T1 = float(sys.argv[1])
	T2 = float(sys.argv[2])
	dT = float(sys.argv[3])
except IndexError, ValueError:
	print("Put T1, T2 and dT")

Np = int(round((T2-T1)/dT))

T_lin = np.linspace(T1, T2, Np +1)
print(T_lin)

for T_ in T_lin:
	Nt = int(10**7)
	T = T_
	Ex = T/20
	seed = np.random.randint(10**9)
	print("./glatzprogram triJumps=0 outpre=../dataMott2/ writelines=10000 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0 " + " T=" + str(T) + " seed2=" + str(seed) + " screen=0.5")
	os.system("./glatzprogram triJumps=0 outpre=../dataMott2/ writelines=10000 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0" + " T=" + str(T) + " seed2=" + str(seed) + " screen=0.5")

