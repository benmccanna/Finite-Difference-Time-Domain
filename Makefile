GCCFLAGS=-Wall -std=c++11

fdtd.out: fdtd.cpp scenario.cpp
	g++ $(GCCFLAGS) -o fdtd.out fdtd.cpp
	
.PHONY: images clean
images:
	 $(MAKE) -C output
	 
clean:
	rm fdtd.out
	rm output/*.png
	rm output/*.dat
	rm output/*.eps