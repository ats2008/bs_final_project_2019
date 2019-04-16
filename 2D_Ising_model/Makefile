INPATH 		= include


CXX 		= g++
CXXFLAGS 	= -Wall -std=c++11 -I$(INPATH)
DEPS 		= main.o isinglib.o iolib.o
DEPS_try	= try.o isinglib.o iolib.o stattools.o
DEPS_check_data_file	= check_data_file.o isinglib.o iolib.o stattools.o
DEPS_resampling	= resampling.o isinglib.o iolib.o stattools.o
DEPS_resample_last	= resample_last.o isinglib.o iolib.o stattools.o
DEPS1		= ising_sweep.o isinglib.o iolib.o
DEPS2		= potts_sweep.o isinglib.o iolib.o
DEPS3		= multipleHistogram.o isinglib.o iolib.o stattools.o
TARGET		= main.exe

all: $(TARGET)

multipleHistogram: $(DEPS3)
	$(CXX) $(CXXFLAGS) -o multipleHistogram.exe $(DEPS3)
	
try: $(DEPS_try)
	$(CXX) $(CXXFLAGS) -o try.exe $(DEPS_try)

check_data_file: $(DEPS_check_data_file)
	$(CXX) $(CXXFLAGS) -o check_data_file.exe $(DEPS_check_data_file)
	
resampling: $(DEPS_resampling)
	$(CXX) $(CXXFLAGS) -o resampling.exe $(DEPS_resampling)
	
resample_last: $(DEPS_resample_last)
	$(CXX) $(CXXFLAGS) -o resample_last.exe $(DEPS_resample_last)



ising_sweep: $(DEPS1)
	$(CXX) $(CXXFLAGS) -o ising_sweep.exe $(DEPS1)
	
potts_sweep: $(DEPS2)
	$(CXX) $(CXXFLAGS) -o potts_sweep.exe $(DEPS2)

$(TARGET) : $(DEPS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(DEPS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
try.o: try.cpp
	$(CXX) $(CXXFLAGS) -c try.cpp
check_data_file.o: check_data_file.cpp
	$(CXX) $(CXXFLAGS) -c check_data_file.cpp
resampling.o: resampling.cpp
	$(CXX) $(CXXFLAGS) -c resampling.cpp
resample_last.o: resample_last.cpp
	$(CXX) $(CXXFLAGS) -c resample_last.cpp
isinglib.o: isinglib.cpp
	$(CXX) $(CXXFLAGS) -c isinglib.cpp
iolib.o: iolib.cpp
	$(CXX) $(CXXFLAGS) -c iolib.cpp
ising_sweep.o: ising_sweep.cpp
	$(CXX) $(CXXFLAGS) -c ising_sweep.cpp
potts_sweep.o: potts_sweep.cpp
	$(CXX) $(CXXFLAGS) -c potts_sweep.cpp
stattools.o: stattools.cpp
	$(CXX) $(CXXFLAGS) -c stattools.cpp
multipleHistogram.o: multipleHistogram.cpp
	$(CXX) $(CXXFLAGS) -c multipleHistogram.cpp
	
clean :
	@rm *.o
	@rm *.exe
	@rm *~
	@rm include/*~
