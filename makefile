all: serial parallel
clean:
	rm serial parallel
serial: serial.cpp
	g++ -o serial serial.cpp
parallel: parallel.cpp
	mpiCC -o parallel parallel.cpp
