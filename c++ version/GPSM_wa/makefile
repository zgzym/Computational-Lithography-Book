main:main.o function.o
	g++ -g -o main main.o function.o -std=c++11 `pkg-config --libs --cflags opencv4`

function.o:function.cc
	g++ -g -c function.cc -std=c++11 `pkg-config --libs --cflags opencv4`

main.o:main.cc
	g++ -g -c main.cc -std=c++11 `pkg-config --libs --cflags opencv4`

.PHONY:clean
clean:
	-rm -f main.o function.o
