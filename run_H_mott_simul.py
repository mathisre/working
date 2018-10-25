import os
import sys
import numpy as np


try:
	T = float(sys.argv[1])
	H1 = float(sys.argv[2])
	H2 = float(sys.argv[3])
	dH = float(sys.argv[4])
except IndexError, ValueError:
	print("Put T, H1 , H2 and dH")
	exit()

Np = int(round((H2-H1)/dH))
print(Np)

H_lin = np.linspace(H1, H2, Np+1)
assert ((H_lin[1]-H_lin[0]) <= dH + 10**-4) and ((H_lin[1]-H_lin[0]) >= dH - 10**-4)
print(H_lin)

Ex = T/10

for Hz in H_lin:
	Nt = int(2*10**7)
	seed = np.random.randint(10**6)
	seedconfig = np.random.randint(9*10**5)
	seedstate = np.random.randint(9*10**5)
	
	print("./glatzprogram outpre=../hallData/ Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=" + str(Hz) + " T=" + str(T) + " seed2=" + str(seed) + " rsconfig=" + str(seedconfig) + " rsstate=" + str(seedstate) + " screen=0.5")
	os.system("./glatzprogram outpre=../hallData/ Nt=" + str(Nt) + " Ex=" + str(Ex) + " Hz=" + str(Hz) + " T=" + str(T) + " seed2=" + str(seed) + " rsconfig=" + str(seedconfig) + " rsstate=" + str(seedstate) + " screen=0.5")

