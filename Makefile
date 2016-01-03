GCCFLAGS=-Wall

fdtd.out: fdtd.cpp
	g++ $(GCCFLAGS) -o fdtd.out fdtd.cpp