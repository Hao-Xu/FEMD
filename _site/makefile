# FORTRAN = gfortran
CC = g++
SOURCES = main.cpp Modeldsid.cpp
DFLAG = -g
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = dsid.x

all : $(EXECUTABLE) $(OBJECTS)

debug: 
	${CC} ${DFLAG} $(SOURCES) -o $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ 

$(OBJECTS):
	$(CC) $(SOURCES) -c   

clean :
	rm -f *.out $(EXECUTABLE) *.o

