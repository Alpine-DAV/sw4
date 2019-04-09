# For MacOS
CXX = mpicxx
CXXFLAGS  += -g -Xpreprocessor -fopenmp -Wno-c++11-extensions 
LINKFLAGS += -framework Accelerate -lomp 

SW4_INCLUDE = /Users/hans/repos/sw4/src
SW4_OBJECTS = /Users/hans/repos/sw4/debug_mp

# For hdf5
HDF5ROOT   = /usr/local/Cellar/hdf5/1.10.2_1
CXXFLAGS  += -I$(HDF5ROOT)/include -DUSE_HDF5 -I$(SW4_INCLUDE)/double -I$(SW4_INCLUDE)
LINKFLAGS += $(SW4_OBJECTS)/*.o -L$(HDF5ROOT)/lib -lhdf5 -lgfortran -lomp -lproj

rfile_2_sfile:
	/bin/rm -f $(SW4_OBJECTS)/main.o;
	$(CXX) $(CXXFLAGS) driver.C \
		-o rfile_2_sfile.exe $(LINKFLAGS) $(LINKLIBS) \

clean:
	/bin/rm -f *.o; /bin/rm -f rfile_2_sfile.exe


