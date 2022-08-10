all:
	g++ -std=c++11 QDMLL_2S.cpp -o QDMLL -Ofast -L $$PWD/incl/ -lfftw3 -lm -Wall
