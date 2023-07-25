all:
	g++ -std=c++11 -O0 -ggdb3 -I../../src -Wall xor.cpp -o xor
	g++ -std=c++11 -O0 -ggdb3 -I../../src -Wall xor_independent.cpp -o xor_independent
	g++ -std=c++11 -DCHANGEABLE_ACTIVATION_AND_AGGREGATION -O0 -ggdb3 -I../../src -Wall xor_changeable_activation_aggregation.cpp -o xor_changeable_activation_aggregation

clean:
	rm xor; rm xor_independent; rm fit\ \=\ *; rm xor_test.res;
