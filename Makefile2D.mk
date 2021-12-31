CC = /usr/local/bin/g++-9
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I/usr/local/opt/icu4c/include -I/usr/local/include#-ffast-math -Ofast -ffinite-math-only
LDFLAGS		=-fopenmp -std=c++17 -L/usr/local/opt/icu4c/lib -I/usr/local/include#-Ofast
SOURCES		=./examples/tutorial.cpp
# SOURCES		=./examples/testDAFMM2D.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testDAFMM2D
export CPATH = /usr/local/include
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
	rm a.out testDAFMM2D *.o
