CXX = mpic++
OPENMPFLAG = -fopenmp
OPENMPLINKINGFLAG =
LIBS = #-Linclude/lis-2.0.32/BUILD/lib -llis
CXXFLAGS = -std=c++17 -g -O2 -fpic 
CPPFLAGS = -DNO_GPU -DNDEBUG -Iinclude/eigen3 -Iinclude -Iinclude/eigen_forma -Iinclude/eigen_forma/Utilities -DENABLE_TIMING 
OBJ = main_final_H.o utils_H.o


.PHONY : all clean distclean

all : main


main : $(OBJ)
	$(CXX) $(OPENMPLINKINGFLAG) $(LIBS) $(OBJ) $(OPENMPFLAG) -o main.exe

$(OBJ) : %.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMPFLAG) -c $<

clean :
	$(RM) *.o *.exe

distclean : clean
	$(RM) main.exe

