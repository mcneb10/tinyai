#ifndef _ARTIFICIAL_NEURAL_NETWORK_HPP_
#define _ARTIFICIAL_NEURAL_NETWORK_HPP_

#include <unordered_map>
#include <cmath>
#include <cstdio>
#include <array>
#include <stack>
#if !defined(INDEPENDENT)
#include <iostream>
#include <fstream>
#endif
#include <algorithm>
#include <vector>

#include "tinyneat.hpp"

namespace ann
{

	enum type
	{
		RECURRENT,
		NON_RECURRENT
	};

	class neuron
	{
	public:
		int type = 0; // 0 = ordinal, 1 = input, 2 = output, 3 = bias
		double value = 0.0;
		bool visited = false;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		double (*aggregation)(std::vector<double>) = neat::aggregation::sum;
		std::string aggregation_name = "sum";
		double (*activation)(double) = neat::activation::sigmoid;
		std::string activation_name = "sigmoid";
#endif
		std::vector<std::pair<size_t, double>> in_nodes;
		neuron() {}
		~neuron() { in_nodes.clear(); }
	};

	// NOTE: you MUST set activation_funcs and aggregation_funcs when you make a ann::neuralnet
	// otherwise things will break!
	class neuralnet
	{
	private:
		std::vector<neuron> nodes;
		bool recurrent = false;

		std::vector<size_t> input_nodes;
		std::vector<size_t> bias_nodes;
		std::vector<size_t> output_nodes;
#if !defined(CHANGEABLE_ACTIVATION_AND_AGGREGATION)
		double sigmoid(double x)
		{
			return 2.0 / (1.0 + std::exp(-4.9 * x)) - 1;
		}
#else
		std::map<std::string, double (*)(double)> activation_funcs;
		std::map<std::string, double (*)(std::vector<double>)> aggregation_funcs;
#endif
		void evaluate_nonrecurrent(const std::vector<double> &input, std::vector<double> &output)
		{

			for (size_t i = 0; i < nodes.size(); i++)
				nodes[i].value = 0.0, nodes[i].visited = false;

			for (size_t i = 0; i < input.size() && i < input_nodes.size(); i++)
			{
				nodes[input_nodes[i]].value = input[i];
				nodes[input_nodes[i]].visited = true;
			}

			for (size_t i = 0; i < bias_nodes.size(); i++)
			{
				nodes[bias_nodes[i]].value = 1.0;
				nodes[bias_nodes[i]].visited = true;
			}

			std::stack<size_t> s;
			for (size_t i = 0; i < output_nodes.size(); i++)
				s.push(output_nodes[i]);

			while (!s.empty())
			{
				size_t t = s.top();

				if (nodes[t].visited == true)
				{
// first is the other node
// second is weight
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
					std::vector<double> input;
					for (size_t i = 0; i < nodes[t].in_nodes.size(); i++)
						input.push_back(nodes[nodes[t].in_nodes[i].first].value * nodes[t].in_nodes[i].second);
					nodes[t].value = nodes[t].activation(nodes[t].aggregation(input));
#else
double sum = 0.0;
					for (size_t i = 0; i < nodes[t].in_nodes.size(); i++)
						sum += nodes[nodes[t].in_nodes[i].first].value * nodes[t].in_nodes[i].second;
					nodes[t].value = sigmoid(sum);
#endif
					s.pop();
				}

				else
				{
					nodes[t].visited = true;

					for (size_t i = 0; i < nodes[t].in_nodes.size(); i++)
						if (nodes[nodes[t].in_nodes[i].first].visited == false)
							// if we haven't calculated value for this node
							s.push(nodes[t].in_nodes[i].first);
				}
			}

			for (size_t i = 0; i < output_nodes.size() && i < output.size(); i++)
				output[i] = nodes[output_nodes[i]].value;
		}
		// FIXME: this probably isn't implemented correctly
		void evaluate_recurrent(const std::vector<double> &input, std::vector<double> &output)
		{

			for (size_t i = 0; i < input.size() && i < input_nodes.size(); i++)
			{
				nodes[input_nodes[i]].value = input[i];
				nodes[input_nodes[i]].visited = true;
			}

			for (size_t i = 0; i < bias_nodes.size(); i++)
			{
				nodes[bias_nodes[i]].value = 1.0;
				nodes[bias_nodes[i]].visited = true;
			}

			// in non-recurrent, each node we will visit only one time per
			// simulation step (similar to the real world)
			// and the values will be saved till the next simulation step
			for (size_t i = 0; i < nodes.size(); i++)
			{
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
	std::vector<double> input;
	for (size_t j = 0; j < nodes[i].in_nodes.size(); j++) //FIXME: is this supposed to be addition
					input.push_back(nodes[nodes[i].in_nodes[j].first].value + nodes[i].in_nodes[j].second);
				if (nodes[i].in_nodes.size() > 0)
					nodes[i].value = nodes[i].activation(nodes[i].aggregation(input));
#else
double sum = 0.0;
				for (size_t j = 0; j < nodes[i].in_nodes.size(); j++)
					sum += nodes[nodes[i].in_nodes[j].first].value + nodes[i].in_nodes[j].second;
				if (nodes[i].in_nodes.size() > 0)
					nodes[i].value = sigmoid(sum);
#endif
			}

			for (size_t i = 0; i < output_nodes.size() && i < output.size(); i++)
				output[i] = nodes[output_nodes[i]].value;
		}

	public:
		neuralnet() {}

		void from_genome(const neat::genome &a)
		{

			unsigned int input_size = a.network_info.input_size;
			unsigned int output_size = a.network_info.output_size;
			unsigned int bias_size = a.network_info.bias_size;

			this->recurrent = a.network_info.recurrent;

			nodes.clear();
			input_nodes.clear();
			bias_nodes.clear();
			output_nodes.clear();

			// We can ignore these nodes as they don't use the ac/ag functions?

			neuron tmp;
			for (unsigned int i = 0; i < input_size; i++)
			{
				nodes.push_back(tmp);
				nodes.back().type = 1;
				this->input_nodes.push_back(nodes.size() - 1);
			}
			for (unsigned int i = 0; i < bias_size; i++)
			{
				nodes.push_back(tmp);
				nodes.back().type = 3;
				this->bias_nodes.push_back(nodes.size() - 1);
			}
			for (unsigned int i = 0; i < output_size; i++)
			{
				nodes.push_back(tmp);
				nodes.back().type = 2;
				this->output_nodes.push_back(nodes.size() - 1);
			}

			std::map<unsigned int, unsigned int> table;
			for (unsigned int i = 0;
				 i < input_nodes.size() + output_nodes.size() + bias_nodes.size(); i++)
				table[i] = i;

			for (auto it = a.genes.begin(); it != a.genes.end(); it++)
			{
				if (!(*it).second.enabled)
					continue;

				neuron n;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
				n.aggregation = it->second.aggregation;
				n.activation = it->second.activation;
				n.aggregation_name = it->second.aggregation_name;
				n.activation_name = it->second.activation_name;
#endif

				if (table.find((*it).second.from_node) == table.end())
				{
					nodes.push_back(n);
					table[(*it).second.from_node] = nodes.size() - 1;
				}
				if (table.find((*it).second.to_node) == table.end())
				{
					nodes.push_back(n);
					table[(*it).second.to_node] = nodes.size() - 1;
				}
			}

			for (auto it = a.genes.begin(); it != a.genes.end(); it++)
				nodes[table[(*it).second.to_node]].in_nodes.push_back(
					std::make_pair(table[(*it).second.from_node], (*it).second.weight));
		}

		void evaluate(const std::vector<double> &input, std::vector<double> &output)
		{
			if (recurrent)
				this->evaluate_recurrent(input, output);
			else
				this->evaluate_nonrecurrent(input, output);
		}
#ifdef INDEPENDENT
		void import_fromfile(std::virtual_file &file)
		{
			file.rewind();
			this->nodes.clear();
			this->input_nodes.clear();
			this->output_nodes.clear();

			std::string rec = file.get_str();
			if (rec == "recurrent")
				this->recurrent = true;
			if (rec == "non-recurrent")
				// there used to be a semi-bug here that files had this as
				// "non-recurrent" but it was checking for "non_recurrent",
				// but non-recurrent is the default anyway
				this->recurrent = false;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
			rec = file.get_str();
			if (rec != "CHANGEABLE_ACTIVATION_AND_AGGREGATION")
			{
				printf("[ERROR] Failed to import a neural network without CHANGEABLE_ACTIVATION_AND_AGGREGATION enabled");
				exit(-1);
			}
#endif

			unsigned int neuron_number;
			file.scanf("%u", &neuron_number);
			this->nodes.resize(neuron_number);

			for (unsigned int i = 0; i < neuron_number; i++)
			{
				unsigned int input_size, type; // 0 = ordinal, 1 = input, 2 = output
				nodes[i].value = 0.0;
				nodes[i].visited = false;

				file.scanf("%u", &type);
				if (type == 1)
					input_nodes.push_back(i);
				if (type == 2)
					output_nodes.push_back(i);
				if (type == 3)
					bias_nodes.push_back(i);

				nodes[i].type = type;

				file.scanf("%u", &input_size);

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
				nodes[i].aggregation_name = file.get_str();
				nodes[i].activation_name = file.get_str();
				nodes[i].aggregation = this->aggregation_funcs.at(nodes[i].aggregation_name);
				nodes[i].activation = this->activation_funcs.at(nodes[i].activation_name);
#endif
				for (unsigned int j = 0; j < input_size; j++)
				{
					unsigned int t; // connected node
					double w;		// weight
					file.scanf("%u", &t);
					file.scanf("%lf", &w);
					nodes[i].in_nodes.push_back(std::make_pair(t, w));
				}
			}
		}
#else
		void import_fromfile(std::string filename)
		{
			std::ifstream o;
			o.open(filename);

			this->nodes.clear();
			this->input_nodes.clear();
			this->output_nodes.clear();

			try
			{
				if (!o.is_open())
					throw "error: cannot open file!";

				std::string rec;
				o >> rec;
				if (rec == "recurrent")
					this->recurrent = true;
				if (rec == "non-recurrent")
					this->recurrent = false;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
				o >> rec;
				if (rec != "CHANGEABLE_ACTIVATION_AND_AGGREGATION")
				{
					throw "Failed to import a neural network without CHANGEABLE_ACTIVATION_AND_AGGREGATION enabled";
				}
#endif

				unsigned int neuron_number;
				o >> neuron_number;
				this->nodes.resize(neuron_number);

				for (unsigned int i = 0; i < neuron_number; i++)
				{
					unsigned int input_size, type; // 0 = ordinal, 1 = input, 2 = output
					nodes[i].value = 0.0;
					nodes[i].visited = false;

					o >> type;
					if (type == 1)
						input_nodes.push_back(i);
					if (type == 2)
						output_nodes.push_back(i);
					if (type == 3)
						bias_nodes.push_back(i);

					nodes[i].type = type;

					o >> input_size;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
					o >> nodes[i].aggregation_name;
					o >> nodes[i].activation_name;
					nodes[i].aggregation = this->aggregation_funcs.at(nodes[i].aggregation_name);
					nodes[i].activation = this->activation_funcs.at(nodes[i].activation_name);
#endif

					for (unsigned int j = 0; j < input_size; j++)
					{
						unsigned int t;
						double w;
						o >> t >> w;
						nodes[i].in_nodes.push_back(std::make_pair(t, w));
					}
				}
			}
			catch (std::string error_message)
			{
				std::cerr << error_message << "\n";
			}

			o.close();
		}
#endif

// NOTE: std::endl replace with \n for more cross platform-ness
#if INDEPENDENT
		void export_tofile(std::virtual_file &file)
		{
			if (this->recurrent)
				file.concat("recurrent\n");
			else
				file.concat("non-recurrent\n");

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
			file.concat("CHANGEABLE_ACTIVATION_AND_AGGREGATION\n");
#endif

			file.concat(std::to_string(nodes.size()) + "\n\n");

			for (size_t i = 0; i < nodes.size(); i++)
			{
				file.concat(std::to_string(nodes[i].type) + " " + std::to_string(nodes[i].in_nodes.size()) + "\n");
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
				file.concat(nodes[i].aggregation_name + "\n");
				file.concat(nodes[i].activation_name + "\n");
#endif
				for (unsigned int j = 0; j < nodes[i].in_nodes.size(); j++)
				{
					file.concat(std::to_string(nodes[i].in_nodes[j].first) + " " +
								std::to_string(nodes[i].in_nodes[j].second) + " ");
				}
				file.concat("\n\n");
			}
		}
#else
		void export_tofile(std::string filename)
		{
			std::ofstream o;
			o.open(filename);

			if (this->recurrent)
				o << "recurrent\n";
			else
				o << "non-recurrent\n";

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
			o << "CHANGEABLE_ACTIVATION_AND_AGGREGATION\n";
#endif

			o << nodes.size() << "\n\n";

			for (size_t i = 0; i < nodes.size(); i++)
			{
				o << nodes[i].type << " ";
				o << nodes[i].in_nodes.size() << "\n";
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
				o << nodes[i].aggregation_name << "\n";
				o << nodes[i].activation_name << "\n";
#endif
				for (unsigned int j = 0; j < nodes[i].in_nodes.size(); j++)
					o << nodes[i].in_nodes[j].first << " "
					  << nodes[i].in_nodes[j].second << " ";
				o << "\n\n";
			}
			o.close();
		}
#endif
	};

} // end of namespace ann

#endif
