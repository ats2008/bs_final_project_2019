INPATH 		= include
CXX 		= g++-8
CXXFLAGS 	= -Wall -fopenmp -std=c++11 -g -lgslcblas -lm -lgsl -I$(INPATH) -I /home/aravind/boost_1_69_0/
DEPS 		= main.o latticelib.o iolib.o
DEPS_resampling	= resampling.o latticelib.o stattools.o
DEPS_resampling_iid	= resampling_idd.o latticelib.o stattools.o
DEPS_stat_analysis   = stat_analysis.o latticelib.o stattools.o 
DEPS_multipleHistogram = multipleHistogram.o latticelib.o stattools.o 
DEPS_joinLatticeHistories	= joinLatticeHistories.o latticelib.o stattools.o
DEPS_joinLatticeProgress	= joinLatticeProgress.o latticelib.o stattools.o
DEPS_omp_thermalize	= omp_thermalize_lattice.o latticelib.o stattools.o
DEPS_omp_progress	= omp_progress_lattice.o latticelib.o stattools.o
TARGET		= main.exe

all: $(TARGET)
	
resampling: $(DEPS_resampling)
	$(CXX) $(CXXFLAGS) -o resampling.exe $(DEPS_resampling)
resampling_idd: $(DEPS_resampling_iid)
	$(CXX) $(CXXFLAGS) -o resampling_idd.exe $(DEPS_resampling_iid)
stat_analysis: $(DEPS_stat_analysis)
	$(CXX) $(CXXFLAGS) -o stat_analysis.exe $(DEPS_stat_analysis)
joinLatticeHistories: $(DEPS_joinLatticeHistories)
	$(CXX) $(CXXFLAGS) -o joinLatticeHistories.exe $(DEPS_joinLatticeHistories)
joinLatticeProgress: $(DEPS_joinLatticeProgress)
	$(CXX) $(CXXFLAGS) -o joinLatticeProgress.exe $(DEPS_joinLatticeProgress)
omp_thermalize_lattice: $(DEPS_omp_thermalize)
	$(CXX) $(CXXFLAGS) -o omp_thermalize_lattice.exe $(DEPS_omp_thermalize)
omp_progress_lattice: $(DEPS_omp_progress)
	$(CXX) $(CXXFLAGS) -o omp_progress_lattice.exe $(DEPS_omp_progress)
multipleHistogram: $(DEPS_multipleHistogram)
	$(CXX) $(CXXFLAGS) -o multipleHistogram.exe $(DEPS_multipleHistogram)
$(TARGET) : $(DEPS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(DEPS)
main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c src/main.cpp
resampling.o: src/resampling.cpp
	$(CXX) $(CXXFLAGS) -c src/resampling.cpp
resampling_idd.o: src/resampling_idd.cpp
	$(CXX) $(CXXFLAGS) -c src/resampling_idd.cpp
joinLatticeHistories.o: src/joinLatticeHistories.cpp
	$(CXX) $(CXXFLAGS) -c src/joinLatticeHistories.cpp
joinLatticeProgress.o: src/joinLatticeProgress.cpp
	$(CXX) $(CXXFLAGS) -c src/joinLatticeProgress.cpp
latticelib.o: src/latticelib.cpp
	$(CXX) $(CXXFLAGS) -c src/latticelib.cpp
omp_thermalize_lattice.o: src/omp_thermalize_lattice.cpp
	$(CXX) $(CXXFLAGS) -c src/omp_thermalize_lattice.cpp
omp_progress_lattice.o: src/omp_progress_lattice.cpp
	$(CXX) $(CXXFLAGS) -c src/omp_progress_lattice.cpp
stattools.o: src/stattools.cpp
	$(CXX) $(CXXFLAGS) -c src/stattools.cpp
stat_analysis.o: src/stat_analysis.cpp
	$(CXX) $(CXXFLAGS) -c src/stat_analysis.cpp
multipleHistogram.o: src/multipleHistogram.cpp
	$(CXX) $(CXXFLAGS) -c src/multipleHistogram.cpp
clean :
	@rm *.o
	@rm *.exe
	@rm *~
