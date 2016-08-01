
# use compiler with openMP support
CC = g++-mp-4.8

# Use pkg-config to lookup the proper compiler and linker flags for LCM
CFLAGS=`pkg-config --cflags lcm` -g -fopenmp -std=gnu++11 -Wno-literal-suffix -Wno-narrowing -O3
LDFLAGS=`pkg-config --libs lcm` -g -fopenmp -std=gnu++11 -Wno-literal-suffix -Wno-narrowing -O3

msg_types_dir = ./types/
msg_types=$(msg_types_dir)*.lcm


all: send-listener-async

send-listener-async: exlcm_location_t.o send-listener-async.o 
	$(CC) -o $@ $^ $(LDFLAGS)

# prevent auto-generated lcm .c/.h files from being deleted
.SECONDARY : exlcm_location_t.c exlcm_location_t.h

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< 

exlcm_%.c exlcm_%.h:
	lcm-gen -c $(msg_types)

clean:
	rm -f send-listener-async
	rm -f *.o
	rm -f exlcm_*
