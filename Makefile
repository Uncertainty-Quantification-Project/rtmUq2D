# $(CC) $(CFLAGS) $(RTMDIR)/*.c  -o $@
# $(CC) $(CFLAGS) $(INCLUDE) $(ZFPDIR)/*.c $(CLIBS) -o $@

CC=mpiicc

RTMDIR = ./SRC
ZFPDIR = ./ZFP_LIBRARY_0.5.5/build/src
TARGETS = rtmuq2d
LIBS = -L./ZFP_LIBRARY_0.5.5/build/lib -lzfp -Wl,-rpath,./ZFP_LIBRARY_0.5.5/build/lib
CLIBS = $(LIBS) -lm
INCLUDE = -I./ZFP_LIBRARY_0.5.5/include
# CFLAGS = -O3 -fopenmp -xHost -inline-forceinline -unroll-aggressive -qopt-assume-safe-padding
CFLAGS = -O3 -fopenmp -xHost

all: $(TARGETS)

rtmuq2d: ./ZFP_LIBRARY_0.5.5/build/lib/$(LIBZFP) $(RTMDIR)/*.c
	$(CC) $(CFLAGS) $(INCLUDE) $(RTMDIR)/*.c $(CLIBS) -o $@

clean:
	rm -f $(TARGETS)