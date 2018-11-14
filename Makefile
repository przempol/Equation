CC=g++
LD=g++
COPTS=`root-config --cflags` -I/usr/local/root/include -g
LDOPTS=`root-config --libs` -g
CFLAGS=-Wall -pedantic -std=c++11 -O3
AUTHOR= Przemyslaw Nowak
#C++ Files
#SOURCES =  main.cpp  Dict.cpp xmath.cpp charts.cpp
SOURCES =  main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
#Dictionary classes
#HEADERS = xmath.h constants.h charts.h linkdef.h
HEADERS = 
EXECUTABLE=program
all: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS) 
	$(LD) -o $@ $^ $(LDOPTS) $(CFLAGS)
	rm -f *.o
#C++ files
.cpp.o:
	$(CC) -o $@ $^ -c $(COPTS)
#Dictionary for ROOT classes
Dict.cpp: $(HEADERS)
	@echo "Generating dictionary ..."
	@rootcint -f  Dict.cpp -c -P -I$ROOTSYS  $(HEADERS)
clean:;		@rm -f $(OBJECTS)  $(EXECUTABLE) *.o *.d Dict.cpp Dict.h



