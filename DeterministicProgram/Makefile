CXX = g++
OPENMPFLAG = 
OPENMPLINKINGFLAG =
CXXFLAGS = -std=c++11 -g -O2 -fpic 
CPPFLAGS = -DNDEBUG -Iinclude
PROGRAM = main_final_eta
CPPFLAGS = -DNDEBUG -Iinclude -DENABLE_TIMING
OBJ = main_final_H.o utils_H.o


#LIB_CHOLMOD = -Linclude/suitesparse/SuiteSparse/CHOLMOD/lib -lcholmod #-Linclude/suitesparse/SuiteSparse/lib
#LIB_AMD = -Linclude/suitesparse/SuiteSparse/AMD/lib -lamd
#LIB_CAMD = -Linclude/suitesparse/SuiteSparse/CAMD/lib -lcamd
#LIB_COLAMD = -Linclude/suitesparse/SuiteSparse/COLAMD/lib -lcolamd
#LIB_CCOLAMD = -Linclude/suitesparse/SuiteSparse/CCOLAMD/lib -lccolamd
#LIB_SUITESPARSECONFIG = -Linclude/suitesparse/SuiteSparse/SuiteSparse_config -lsuitesparseconfig



.PHONY : all clean distclean

all : main


main : $(OBJ)
	$(CXX) $(OPENMPLINKINGFLAG) $(OBJ) -o main.exe


$(OBJ) : %.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMPFLAG) -c $<



clean :
	$(RM) *.o

distclean : clean
	$(RM) main.exe

