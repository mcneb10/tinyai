#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <tinyneat.hpp>
#include <tinyann.hpp>

// returns the fitness.
unsigned int xor_test(ann::neuralnet& n, bool write_output){
	std::vector<double> input(2, 0.0);
	std::vector<double> output(1, 0.0);
	unsigned int fitness = 0;
	double answer;

	if (write_output) std::cerr << "     > begin xor test" << std::endl << "        (";

	input[0] = 0.0, input[1] = 0.0, answer = 0.0;
	n.evaluate(input, output);
	if (write_output) std::cerr << output[0] << " ";
	fitness += std::min(1.0 / ((answer - output[0]) * (answer - output[0])), 50.0);	

	input[0] = 0.0, input[1] = 1.0, answer = 1.0;
	n.evaluate(input, output);
	if (write_output) std::cerr << output[0] << " ";
	fitness += std::min(1.0 / ((answer - output[0]) * (answer - output[0])), 50.0);	

	input[0] = 1.0, input[1] = 0.0, answer = 1.0;
	n.evaluate(input, output);
	if (write_output) std::cerr << output[0] << " ";
	fitness += std::min(1.0 / ((answer - output[0]) * (answer - output[0])), 50.0);	

	input[0] = 1.0, input[1] = 1.0, answer = 0.0;
	n.evaluate(input, output);
	if (write_output) std::cerr << output[0] << ")";
	fitness += std::min(1.0 / ((answer - output[0]) * (answer - output[0])), 50.0);	

	if (write_output) std::cerr << " " << fitness << std::endl;

	return fitness;
}


void test_output(std::vector<std::virtual_file> files){
	for (std::virtual_file file : files)
	{
		if (file.name == "fitc = 200")
		{
			ann::neuralnet n;
			n.activation_funcs = {{std::string("sigmoid"), &neat::activation::sigmoid}};
			n.aggregation_funcs = {{std::string("sum"), &neat::aggregation::sum}};
			n.import_fromfile(file);
			xor_test(n, true);
			return;
		}
	}
}

int main(){
	std::vector<std::virtual_file> files;
	neat::pool p(2, 1, 0, false);
	//p.import_fromfile("xor_test.resc");	
	srand(time(NULL));
	unsigned int max_fitness = 0;
	while (max_fitness < 200){
		unsigned int current_fitness = 0;
		unsigned int min_fitness = 100000;	
		for (auto s = p.species.begin(); s != p.species.end(); s++)
			for (size_t i=0; i<(*s).genomes.size(); i++){				
				ann::neuralnet n;
				neat::genome& g = (*s).genomes[i];
				n.from_genome(g);
				current_fitness = xor_test(n, false);
				if (current_fitness < min_fitness)
					min_fitness = current_fitness;
				if (current_fitness > max_fitness){
					max_fitness = current_fitness;
					std::virtual_file file("fitc = " + std::to_string(current_fitness));
					n.export_tofile(file);
					files.push_back(file);
				}
				g.fitness = current_fitness;
			}

		std::cerr << "Generation " << p.generation() << " successfuly tested. Global min fitness: " << min_fitness << ", Global max fitness: " << max_fitness << std::endl;
		p.new_generation();
	}

	test_output(files);
	std::virtual_file result("xor_test.resc");
	p.export_tofile(result);
	return 0;
}
