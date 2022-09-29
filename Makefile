# CC: compiler
# LBoost: links to Boost library
# LGSL: links to GSL library
#
# With Boost library and small files:
#
# LBoost=-I ~/Boost/boost_1_57_0
# bosons : main.cpp bosonarray.hpp boson_eom.hpp boson_IC.hpp boson_observables.hpp
# $(CC) $(LBoost) $(LGSL) $(LFFTW) -o bosons main.cpp

CC=g++ -O3 -fopenmp
LGSL=-L/usr/local/lib -I/usr/local/include -lgsl
LFFTW=-lfftw3 -lfftw3_threads -lm -lpthread


bosons_F : main_F.cpp boson_gpe_F.cpp boson_gpe_F.hpp
	$(CC) $(LGSL) $(LFFTW) -o bosons_F main_F.cpp boson_gpe_F.cpp

clean :
	rm bosons_F
