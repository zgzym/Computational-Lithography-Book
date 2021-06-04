OPC_tv:OPC_tv.o function.o
	g++ -g -o OPC_tv OPC_tv.o function.o -std=c++11 `pkg-config --libs --cflags opencv4`

function.o:function.cc
	g++ -g -c function.cc -std=c++11 `pkg-config --libs --cflags opencv4`

OPC_tv.o:OPC_tv.cc
	g++ -g -c OPC_tv.cc -std=c++11 `pkg-config --libs --cflags opencv4`

.PHONY:clean
clean:
	-rm -f OPC_tv.o function.o
