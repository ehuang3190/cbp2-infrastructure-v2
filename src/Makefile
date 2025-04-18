CXX		=	g++
CXXFLAGS	=	-g -O3 -Wall

all:		predict

predict:	predict.cc trace.cc tage_predictor.cc predictor.h branch.h trace.h my_predictor.h tage_predictor.h
		$(CXX) $(CXXFLAGS) -o predict predict.cc trace.cc tage_predictor.cc

clean:
		rm -f predict
