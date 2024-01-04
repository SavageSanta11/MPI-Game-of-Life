all: life-blocking life-nonblocking

CXX = mpicxx
CXXFLAGS = -std=c++11 -O2

life-blocking: life-blocking.C
	$(CXX) $(CXXFLAGS) -o life-blocking $<

life-nonblocking: life-nonblocking.C
	$(CXX) $(CXXFLAGS) -o life-nonblocking $<

clean:
	rm -f life-blocking.o life-blocking  life-nonblocking life-nonblocking.o