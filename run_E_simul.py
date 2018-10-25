import os
import sys
import numpy as np


try:
	T = float(sys.argv[1])
	E1 = float(sys.argv[2])
	E2 = float(sys.argv[3])
	dE = float(sys.argv[4])
except IndexError, ValueError:
	print("Put T, E1 , E2 and dE")
	exit()

Np = int(round((E2-E1)/dE))
print(Np)

E_lin = np.linspace(E1, E2, Np+1)
assert ((E_lin[1]-E_lin[0]) <= dE + 10**-4) and ((E_lin[1]-E_lin[0]) >= dE - 10**-4)
print(E_lin)


for Ex in E_lin:
	Nt = int(10**7)
	seed = np.random.randint(10**9)
	print("./glatzprogram triJumps=0 outpre=../data_const_T/ writelines=10000 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0 " + " T=" + str(T) + " seed2=" + str(seed))
	os.system("./glatzprogram triJumps=0 outpre=../data_const_T/ writelines=10000 Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=0" + " T=" + str(T) + " seed2=" + str(seed))

