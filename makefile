cc := gcc
cpp := g++
cflags = -g -Wfatal-errors -Ilib/lunar/ -Ilib/sat_code/ -Ilib/maidenhead/ -Idata/tles
cppflags = -g -Wfatal-errors -Ilib/lunar/ -Ilib/sat_code/ -Ilib/maidenhead/ -Idata/tles
ldflags = -lm 

lunar_src = lib/lunar/alt_az.cpp lib/lunar/precess.cpp lib/lunar/miscell.cpp lib/lunar/obliquit.cpp lib/lunar/nutation.cpp lib/lunar/date.cpp lib/lunar/snprintf.cpp
lunar_objs = $(patsubst %.cpp, %.o, $(lunar_src))
sat_code_src = lib/sat_code/sgp.cpp lib/sat_code/sgp4.cpp lib/sat_code/sgp8.cpp lib/sat_code/sdp4.cpp lib/sat_code/sdp8.cpp lib/sat_code/deep.cpp lib/sat_code/basics.cpp lib/sat_code/get_el.cpp lib/sat_code/common.cpp lib/sat_code/tle_out.cpp lib/sat_code/observe.cpp
sat_code_objs = $(patsubst %.cpp, %.o, $(sat_code_src))
maidenhead_src = lib/maidenhead/maidenhead.c
maidenhead_objs = $(patsubst %.c, %.o, $(maidenhead_src))

.PHONY: run
run: tarsat
	-rm *.csv
	./tarsat
	gnuplot graphs.gnuplot  -p
	grep '' 25544_*_N.csv

tarsat: main.c satellite.o  ${sat_code_objs} ${lunar_objs}  ${maidenhead_objs}
	${cc} ${cflags} -o $@ $^ ${ldflags}

	#data/tles/generated/sat_tles.c
data/tles/generated/sat_tles.c: data/tles/nasabare.txt data/tles/process.py
	make -C data/tles update
	make -C data/tles run

%.o: %.cpp
	${cpp} ${cppflags} -c -o $@ $^ 
%.o: %.c
	${cc} ${cflags} -c -o $@ $^ 

.PHONY: lunar
lunar: ${lunar_objs}
.PHONY: sat_code
sat_code: ${sat_code_objs}

clean:
	-rm *.csv tarsat satellite.o ${sat_code_objs} ${lunar_objs}
