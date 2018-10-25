CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11
glatzprogram: paramfile.o stringutils.o systemclass.o main.o fileutils.o CGanalysis.o MKcurrent.o params.o run.o treeutils.o mersenne.o
	g++  -Wall  -O2 -o glatzprogram paramfile.o stringutils.o systemclass.o main.o fileutils.o CGanalysis.o MKcurrent.o params.o run.o treeutils.o mersenne.o
paramfile.o:
stringutils.o:
systemclass.o:
main.o:
fileutils.o:
CGanalysis.o:
MKcurrent.o:
params.o:
run.o:
treeutils.o:
mersenne.o:

clean:
	rm *.o
