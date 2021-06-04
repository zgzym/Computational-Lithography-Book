GPSM_tv:GPSM_tv.o function.o
	g++ -g -o GPSM_tv GPSM_tv.o function.o -std=c++11 `pkg-config --libs --cflags opencv4`

function.o:function.cc
	g++ -g -c function.cc -std=c++11 `pkg-config --libs --cflags opencv4`

GPSM_tv.o:GPSM_tv.cc
	g++ -g -c GPSM_tv.cc -std=c++11 `pkg-config --libs --cflags opencv4`

.PHONY:clean
clean:
	-rm -f GPSM_tv.o function.o
