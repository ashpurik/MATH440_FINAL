all: serial
clean:
	rm serial
serial: serial.cpp
	g++ -o serial serial.cpp
