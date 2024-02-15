CXX	:= c++
CXXFLAGS := -Wall -Wextra -MMD -O2 -m64 -g 
INCLUDES := /home/gurobi/10.0.0/linux64/include/  
LIBSINCLUDES := /home/gurobi/10.0.0/linux64/lib/ 
LIBS = gurobi_c++ gurobi100 m

CPPFLAGS += $(foreach includedir,$(INCLUDES),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LIBSINCLUDES),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

OBJECTS := main_command.o tests/updcorrelation.o data/instance.o data/saainstance.o data/allsceinstance.o input/input.o lshaped/mainproblems/abstractmainproblem.o lshaped/subproblems/linearizationsubproblem.o lshaped/lshaped.o solver/saa/saaproblem.o solver/baseline/linearizationproblem.o solver/baseline/expectedproblem.o solver/saa/saatestproblem.o lshaped/mainproblems/linearizationmainproblem.o solver/status.o

DIRS = $(OBJECTS:.o=.d)

EXE:=exe

$(EXE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS) 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	$(RM) *.o *~ $(EXE)
	rm -rf $(OBJECTS)
	rm -rf $(DIRS)
