CXX = mpic++
OPENMPFLAG = -fopenmp
OPENMPLINKINGFLAG =
LIBS = #-Linclude/lis-2.0.32/BUILD/lib -llis
CXXFLAGS = -std=c++11 -g -O2 -fpic 
CPPFLAGS = -DNDEBUG 
PROGRAM = main_final_eta
CPPFLAGS = -DNDEBUG -Iinclude/eigen3 -Iinclude -Iinclude/eigen_forma -Iinclude/eigen_forma/Utilities -DENABLE_TIMING #-Iinclude/lis-2.0.32/BUILD/include
OBJ = main_final_H.o utils_H.o





.PHONY : all clean distclean

all : main


main : $(OBJ)
	$(CXX) $(OPENMPLINKINGFLAG) $(LIBS) $(OBJ) $(OPENMPFLAG) -o main.exe


$(OBJ) : %.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMPFLAG) -c $<



clean :
	$(RM) *.o

distclean : clean
	$(RM) main.exe

