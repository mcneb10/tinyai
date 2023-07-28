#ifndef __TINYNEAT_HPP__
#define __TINYNEAT_HPP__

/* custom defines:
 * INCLUDE_ENABLED_GENES_IF_POSSIBLE  - if during experiment you found that too
 * many genes are disabled, you can use this option. ALLOW_RECURRENCY_IN_NETWORK
 * - allowing recurrent links
 *
 * GIVING_NAMES_FOR_SPECIES           - giving species unique names (need a
 * dictionary with names in a file "specie_names.dict"
 */

#if !defined(INDEPENDENT)
#include <fstream>
#include <iostream>
#endif
#include "virtual_file.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <map>
#include <queue>
#include <random>
#include <string>
#include <vector>

namespace neat
{
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
	namespace aggregation
	{
		double sum(std::vector<double> &a)
		{
			double r = 0;
			for (auto n : a)
				r += n;
			return r;
		}
		double product(std::vector<double> &a)
		{
			double r = 1;
			for (auto n : a)
				r *= n;
			return r;
		}
		double max(std::vector<double> &a) { return *std::max_element(a.begin(), a.end()); }
		double min(std::vector<double> &a) { return *std::min_element(a.begin(), a.end()); }
		double median(std::vector<double> &a)
		{
			if (a.size() % 2 == 0)
			{
				return (a[a.size() / 2] + a[(a.size() / 2) - 1]) / 2;
			}
			return a[a.size() / 2];
		}
		double mean(std::vector<double> &a)
		{
			double r = 0;
			for (auto n : a)
				r += n;
			return r / a.size();
		}
	} // namespace aggregation
	namespace activation
	{
		double sigmoid(double x) { return 2.0 / (1.0 + std::exp(-4.9 * x)) - 1; }
		double tanh(double a)
		{
			return std::tanh(a);
		}
		double sin(double a)
		{
			return std::sin(a);
		}
		// idk, let the user add the rest
	} // namespace activation
#endif
	typedef struct
	{
		double connection_mutate_chance = 0.25;
		double perturb_chance = 0.90;
		double crossover_chance = 0.75;
		double link_mutation_chance = 2.0;
		double node_mutation_chance = 0.50;
		double bias_mutation_chance = 0.40;
		double step_size = 0.1;
		double disable_mutation_chance = 0.4;
		double enable_mutation_chance = 0.2;
#if CHANGEABLE_ACTIVATION_AND_AGGREGATION
		double activation_mutate_chance = 0.0;
		double aggregation_mutate_chance = 0.0;
#endif
#if INDEPENDENT
		void read(std::virtual_file &file);
		void write(std::virtual_file &file, std::string prefix);
#else
		void read(std::ifstream &o);
		void write(std::ofstream &o, std::string prefix);
#endif
	} mutation_rate_container;

	typedef struct
	{
		unsigned int population = 240;
		double delta_disjoint = 2.0;
		double delta_weights = 0.4;
		double delta_threshold = 1.3;
		unsigned int stale_species = 15;

#if INDEPENDENT
		void read(std::virtual_file &file);
		void write(std::virtual_file &file, std::string prefix);
#else
		void read(std::ifstream &o);
		void write(std::ofstream &o, std::string prefix);
#endif
	} speciating_parameter_container;

	typedef struct
	{
		unsigned int input_size;
		unsigned int bias_size;
		unsigned int output_size;
		unsigned int functional_nodes;
		bool recurrent;
	} network_info_container;

	typedef struct
	{
		unsigned int innovation_num = -1;
		unsigned int from_node = -1;
		unsigned int to_node = -1;
		double weight = 0.0;
		bool enabled = true;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		double (*aggregation)(std::vector<double> &);
		std::string aggregation_name;
		double (*activation)(double);
		std::string activation_name;
#endif
	} gene;

	class genome
	{
	private:
		genome(){};

	public:
		unsigned int fitness = 0;
		unsigned int adjusted_fitness = 0;
		unsigned int global_rank = 0;
		unsigned int max_neuron;
		unsigned int can_be_recurrent = false;

		mutation_rate_container mutation_rates;
		network_info_container network_info;

		std::map<unsigned int, gene> genes;

		genome(network_info_container &info, mutation_rate_container &rates)
		{
			mutation_rates = rates;
			network_info = info;
			max_neuron = network_info.functional_nodes;
		}

		genome(const genome &) = default;
	};

	/* a specie is group of genomes which differences is smaller than some threshold
	 */
	typedef struct
	{
		unsigned int top_fitness = 0;
		unsigned int average_fitness = 0;
		unsigned int staleness = 0;

#ifdef GIVING_NAMES_FOR_SPECIES
		std::string name;
#endif
		std::vector<genome> genomes;
	} specie;

	class innovation_container
	{
	private:
		unsigned int _number;
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> track;
		void set_innovation_number(unsigned int num)
		{
			_number = num;
			reset();
		}
		friend class pool;

	public:
		innovation_container() : _number(0) {}
		void reset() { track.clear(); };
		unsigned int add_gene(gene &g)
		{
			auto it = track.find(std::make_pair(g.from_node, g.to_node));
			if (it == track.end())
				return track[std::make_pair(g.from_node, g.to_node)] = ++_number;
			else
				return (*it).second;
		}
		unsigned int number() { return _number; }
	};

	/* a small world, where individuals (genomes) are making babies and evolving,
	 * becoming better and better after each generation :)
	 */
	class pool
	{
	private:
		pool(){};

		/* important part, only accecible for friend */
		innovation_container innovation;

		/* innovation tracking in current generation, should be cleared after each
		 * generation */
		std::map<std::pair<unsigned int, unsigned int>, unsigned int> track;

		unsigned int generation_number = 1;

		/* evolutionary methods */
		genome crossover(const genome &g1, const genome &g2);
		void mutate_weight(genome &g);
		void mutate_enable_disable(genome &g, bool enable);
		void mutate_link(genome &g, bool force_bias);
		void mutate_node(genome &g);
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		void mutate_activation_func(genome &g);
		void mutate_aggregation_func(genome &g);
#endif
		void mutate(genome &g);

		double disjoint(const genome &g1, const genome &g2);
		double weights(const genome &g1, const genome &g2);
		bool is_same_species(const genome &g1, const genome &g2);

		/* specie ranking */
		void rank_globally();
		void calculate_average_fitness(specie &s);
		unsigned int total_average_fitness();

		/* evolution */
		void cull_species(bool cut_to_one);
		genome breed_child(specie &s);
		void remove_stale_species();
		void remove_weak_species();
		void add_to_species(genome &child);

	public:
		/* pool parameters */
		unsigned int max_fitness = 0;

		/* mutation parameters */
		mutation_rate_container mutation_rates;

		/* species parameters */
		speciating_parameter_container speciating_parameters;

		/* neural network parameters */
		network_info_container network_info;

		// pool's local random number generator
		std::random_device rd;
		std::mt19937 generator;

		/* species */
		std::list<specie> species;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		std::map<std::string, double (*)(double)> activation_funcs;
		std::map<std::string, double (*)(std::vector<double> &)> aggregation_funcs;
		std::string dac;
		std::string dag;
		pool(unsigned int input, unsigned int output, unsigned int bias = 1,
			 bool rec = false, mutation_rate_container m = mutation_rate_container(),
			 speciating_parameter_container s = speciating_parameter_container(),
			 network_info_container c = network_info_container(),
			 std::map<std::string, double (*)(double)> ac = {{std::string("sigmoid"), &neat::activation::sigmoid}},
			 std::map<std::string, double (*)(std::vector<double> &)> ag = {{std::string("sum"), &neat::aggregation::sum}},
			 std::string defaultActivation = "sigmoid", // set these to empty if you want them randomly chosen
			 std::string defaultAggregation = "sum")
		// If default function name null all nodes will have a random aggregation/activation function
		{
#else

		// constructor
		pool(unsigned int input, unsigned int output, unsigned int bias = 1,
			 bool rec = false, mutation_rate_container m = mutation_rate_container(),
			 speciating_parameter_container s = speciating_parameter_container())
		{
#endif
			this->network_info.input_size = input;
			this->network_info.output_size = output;
			this->network_info.bias_size = bias;
			this->network_info.functional_nodes = input + output + bias;
			this->network_info.recurrent = rec;

			this->speciating_parameters = s;
			this->mutation_rates = m;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
			this->activation_funcs = ac;
			this->aggregation_funcs = ag;
			this->dac = defaultActivation;
			this->dag = defaultAggregation;
#endif

			// seed the mersenne twister with
			// a random number from our computer
			generator.seed(rd());

			// create a basic generation with default genomes
			for (unsigned int i = 0; i < this->speciating_parameters.population; i++)
			{
				genome new_genome(this->network_info, this->mutation_rates);
				this->mutate(new_genome);
				this->add_to_species(new_genome);
			}
		}

		/* next generation */
		void new_generation();
		unsigned int generation() { return this->generation_number; }

		/* calculate fitness */
		std::vector<std::pair<specie *, genome *>> get_genomes()
		{
			std::vector<std::pair<specie *, genome *>> genomes;
			for (auto s = this->species.begin(); s != this->species.end(); s++)
				for (size_t i = 0; i < (*s).genomes.size(); i++)
					genomes.push_back(std::make_pair(&(*s), &((*s).genomes[i])));
			return genomes;
		}

/* import and export */
#if INDEPENDENT
		void import_fromfile(std::virtual_file &file);
		void export_tofile(std::virtual_file &file);
#else
		void import_fromfile(std::string filename);
		void export_tofile(std::string filename);
#endif
	};

	/* now the evolutionary functions itself */
	genome pool::crossover(const genome &g1, const genome &g2)
	{
		// Make sure g1 has the higher fitness, so we will include only
		// disjoint/excess genes from the first genome.
		if (g2.fitness > g1.fitness)
			return crossover(g2, g1);
		genome child(this->network_info, this->mutation_rates);

		auto it1 = g1.genes.begin();

		// coin flip random number distributor
		std::uniform_int_distribution<int> coin_flip(1, 2);

		for (; it1 != g1.genes.end(); it1++)
		{

			// if innovation marks match, do the crossover, else include from the first
			// genome because its fitness is not smaller than the second's

			auto it2 = g2.genes.find((*it1).second.innovation_num);
			if (it2 != g2.genes.end())
			{

				// do the coin flip
				int coin = coin_flip(this->generator);

// now, after flipping the coin, we do the crossover.
#ifdef INCLUDE_ENABLED_GENES_IF_POSSIBLE
				if (coin == 2 && (*it2).enabled)
					child.genes[(*it1).second.innovation_num] = (*it2).second;
				else
					child.genes[(*it1).second.innovation_num] = (*it1).second;
#else
				if (coin == 2)
					child.genes[(*it1).second.innovation_num] = (*it2).second;
				else
					child.genes[(*it1).second.innovation_num] = (*it1).second;

#endif
			}
			else
				// as said before, we include the disjoint gene
				// from the first (with larger fitness) otherwise
				child.genes[(*it1).second.innovation_num] = (*it1).second;
		}

		child.max_neuron = std::max(g1.max_neuron, g2.max_neuron);
		return child;
	}

	/* mutations */
	void pool::mutate_weight(genome &g)
	{
		double step = this->mutation_rates.step_size;
		std::uniform_real_distribution<double> real_distributor(0.0, 1.0);

		for (auto it = g.genes.begin(); it != g.genes.end(); it++)
		{
			if (real_distributor(this->generator) < this->mutation_rates.perturb_chance)
				(*it).second.weight +=
					real_distributor(this->generator) * step * 2.0 - step;
			else
				(*it).second.weight = real_distributor(this->generator) * 4.0 - 2.0;
		}
	}

	void pool::mutate_enable_disable(genome &g, bool enable)
	{

		// that shit is safe because there's no changings in map
		// during this function
		std::vector<gene *> v;

		for (auto it = g.genes.begin(); it != g.genes.end(); it++)
			if ((*it).second.enabled != enable)
				v.push_back(&((*it).second));

		if (v.size() == 0)
			return;

		std::uniform_int_distribution<int> distributor(0, v.size() - 1);
		v[distributor(this->generator)]->enabled = enable;
	}

	void pool::mutate_link(genome &g, bool force_bias)
	{
		/* network encoding:
		 * | input nodes | bias | output nodes |
		 */
		auto is_input = [&](unsigned int node) -> bool
		{
			return node < this->network_info.input_size;
		};
		auto is_output = [&](unsigned int node) -> bool
		{
			return node < this->network_info.functional_nodes &&
				   node >=
					   (this->network_info.input_size + this->network_info.bias_size);
		};
		auto is_bias = [&](unsigned int node) -> bool
		{
			return node <
					   (this->network_info.input_size + this->network_info.bias_size) &&
				   node >= this->network_info.input_size;
		};

		std::uniform_int_distribution<unsigned int> distributor1(0, g.max_neuron - 1);
		unsigned int neuron1 = distributor1(this->generator);

		std::uniform_int_distribution<unsigned int> distributor2(
			this->network_info.input_size + this->network_info.bias_size,
			g.max_neuron - 1);
		unsigned int neuron2 = distributor2(this->generator);

		if (is_output(neuron1) && is_output(neuron2))
			return;
		if (is_bias(neuron2))
			return;
		if (neuron1 == neuron2 && (!force_bias))
			return;
		if (is_output(neuron1))
			std::swap(neuron1, neuron2);

		if (force_bias)
		{
			std::uniform_int_distribution<unsigned int> bias_choose(
				this->network_info.input_size,
				this->network_info.input_size + this->network_info.output_size - 1);
			neuron1 = bias_choose(this->generator);
		}

		if (!g.network_info.recurrent)
		{
			// check for recurrency using BFS
			bool has_recurrence = false;
			if (is_bias(neuron1) || is_input(neuron1))
				has_recurrence = false;
			else
			{

				std::queue<unsigned int> que;
				std::vector<std::vector<unsigned int>> connections(g.max_neuron);
				for (auto it = g.genes.begin(); it != g.genes.end(); it++)
					connections[(*it).second.from_node].push_back((*it).second.to_node);
				connections[neuron1].push_back(neuron2);

				for (size_t i = 0; i < connections[neuron1].size(); i++)
					que.push(connections[neuron1][i]);
				while (!que.empty())
				{
					unsigned int tmp = que.front();
					if (tmp == neuron1)
					{
						has_recurrence = true;
						break;
					}
					que.pop();
					for (size_t i = 0; i < connections[tmp].size(); i++)
						que.push(connections[tmp][i]);
				}
			}
			if (has_recurrence)
				return;

			// now we are sure that it doesn't has any recurrency
		}

		// now we can create a link
		gene new_gene;
		new_gene.from_node = neuron1;
		new_gene.to_node = neuron2;

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		if (this->dag.empty())
		{
			auto t = this->aggregation_funcs.begin();
			std::advance(t, rand() % this->aggregation_funcs.size());
			new_gene.aggregation_name = t->first;
		}
		else
		{
			new_gene.aggregation_name = this->dag;
		}
		if (this->dac.empty())
		{
			auto t = this->activation_funcs.begin();
			std::advance(t, rand() % this->activation_funcs.size());
			new_gene.activation_name = t->first;
		}
		else
		{
			new_gene.activation_name = this->dac;
		}
		new_gene.aggregation = this->aggregation_funcs.at(new_gene.aggregation_name);
		new_gene.activation = this->activation_funcs.at(new_gene.activation_name);
#endif

		// if genome already has this connection
		for (auto it = g.genes.begin(); it != g.genes.end(); it++)
			if ((*it).second.from_node == neuron1 && (*it).second.to_node == neuron2)
				return;

		// add new innovation if needed
		new_gene.innovation_num = this->innovation.add_gene(new_gene);

		// mutate new link
		std::uniform_real_distribution<double> weight_generator(0.0, 1.0);
		new_gene.weight = weight_generator(this->generator) * 4.0 - 2.0;

		g.genes[new_gene.innovation_num] = new_gene;
	}

	void pool::mutate_node(genome &g)
	{
		if (g.genes.size() == 0)
			return;

		g.max_neuron++;

		// randomly choose a gene to mutate
		std::uniform_int_distribution<unsigned int> distributor(0,
																g.genes.size() - 1);
		unsigned int gene_id = distributor(this->generator);
		auto it = g.genes.begin();
		std::advance(it, gene_id);

		if ((*it).second.enabled == false)
			return;

		(*it).second.enabled = false;

		gene new_gene1;
		new_gene1.from_node = (*it).second.from_node;
		new_gene1.to_node = g.max_neuron - 1; // to the last created neuron
		new_gene1.weight = 1.0;
		new_gene1.innovation_num = this->innovation.add_gene(new_gene1);
		new_gene1.enabled = true;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		if (this->dag.empty())
		{
			auto t = this->aggregation_funcs.begin();
			std::advance(t, rand() % this->aggregation_funcs.size());
			new_gene1.aggregation_name = t->first;
		}
		else
		{
			new_gene1.aggregation_name = this->dag;
		}
		if (this->dac.empty())
		{
			auto t = this->activation_funcs.begin();
			std::advance(t, rand() % this->activation_funcs.size());
			new_gene1.activation_name = t->first;
		}
		else
		{
			new_gene1.activation_name = this->dac;
		}
		new_gene1.aggregation = this->aggregation_funcs.at(new_gene1.aggregation_name);
		new_gene1.activation = this->activation_funcs.at(new_gene1.activation_name);
#endif

		gene new_gene2;
		new_gene2.from_node = g.max_neuron - 1; // from the last created neuron
		new_gene2.to_node = (*it).second.to_node;
		new_gene2.weight = (*it).second.weight;
		new_gene2.innovation_num = this->innovation.add_gene(new_gene2);
		new_gene2.enabled = true;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		if (this->dag.empty())
		{
			auto t = this->aggregation_funcs.begin();
			std::advance(t, rand() % this->aggregation_funcs.size());
			new_gene2.aggregation_name = t->first;
		}
		else
		{
			new_gene2.aggregation_name = this->dag;
		}
		if (this->dac.empty())
		{
			auto t = this->activation_funcs.begin();
			std::advance(t, rand() % this->activation_funcs.size());
			new_gene2.activation_name = t->first;
		}
		else
		{
			new_gene2.activation_name = this->dac;
		}
		new_gene2.aggregation = this->aggregation_funcs.at(new_gene2.aggregation_name);
		new_gene2.activation = this->activation_funcs.at(new_gene2.activation_name);
#endif

		g.genes[new_gene1.innovation_num] = new_gene1;
		g.genes[new_gene2.innovation_num] = new_gene2;
	}

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
	void pool::mutate_activation_func(genome &g)
	{
		if (g.genes.size() == 0)
			return;

		// randomly choose a gene to mutate
		std::uniform_int_distribution<unsigned int> distributor(0,
																g.genes.size() - 1);
		unsigned int gene_id = distributor(this->generator);
		auto it = g.genes.begin();
		std::advance(it, gene_id);

		auto t = this->activation_funcs.begin();
		std::advance(t, rand() % this->activation_funcs.size());
		it->second.activation_name = it->second.activation_name;
		it->second.activation = this->activation_funcs.at(it->second.activation_name);
	}
	void pool::mutate_aggregation_func(genome &g)
	{
		if (g.genes.size() == 0)
			return;

		// randomly choose a gene to mutate
		std::uniform_int_distribution<unsigned int> distributor(0,
																g.genes.size() - 1);
		unsigned int gene_id = distributor(this->generator);
		auto it = g.genes.begin();
		std::advance(it, gene_id);

		auto t = this->aggregation_funcs.begin();
		std::advance(t, rand() % this->aggregation_funcs.size());
		it->second.aggregation_name = it->second.aggregation_name;
		it->second.aggregation = this->aggregation_funcs.at(it->second.aggregation_name);
	}
#endif

	void pool::mutate(genome &g)
	{
		double coefficient[2] = {0.95, 1.05263};

		std::uniform_int_distribution<int> coin_flip(0, 1);

		g.mutation_rates.enable_mutation_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.disable_mutation_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.connection_mutate_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.node_mutation_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.link_mutation_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.bias_mutation_chance *=
			coefficient[coin_flip(this->generator)];
		g.mutation_rates.crossover_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.perturb_chance *= coefficient[coin_flip(this->generator)];

#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		g.mutation_rates.activation_mutate_chance *= coefficient[coin_flip(this->generator)];
		g.mutation_rates.aggregation_mutate_chance *= coefficient[coin_flip(this->generator)];
#endif

		std::uniform_real_distribution<double> mutate_or_not_mutate(0.0, 1.0);

		if (mutate_or_not_mutate(this->generator) <
			g.mutation_rates.connection_mutate_chance)
			this->mutate_weight(g);

		double p;

		p = g.mutation_rates.link_mutation_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_link(g, false);
			p = p - 1.0;
		}

		p = g.mutation_rates.bias_mutation_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_link(g, true);
			p = p - 1.0;
		}

		p = g.mutation_rates.node_mutation_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_node(g);
			p = p - 1.0;
		}

		p = g.mutation_rates.enable_mutation_chance;
		;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_enable_disable(g, true);
			p = p - 1.0;
		}

		p = g.mutation_rates.disable_mutation_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
				this->mutate_enable_disable(g, false);
			p = p - 1.0;
		}
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		p = this->mutation_rates.aggregation_mutate_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
			{
				this->mutate_aggregation_func(g);
			}
			p = p - 1.0;
		}
		p = this->mutation_rates.activation_mutate_chance;
		while (p > 0.0)
		{
			if (mutate_or_not_mutate(this->generator) < p)
			{
				this->mutate_activation_func(g);
			}
			p = p - 1.0;
		}
#endif
	}

	double pool::disjoint(const genome &g1, const genome &g2)
	{
		auto it1 = g1.genes.begin();
		auto it2 = g2.genes.begin();

		unsigned int disjoint_count = 0;
		for (; it1 != g1.genes.end(); it1++)
			if (g2.genes.find((*it1).second.innovation_num) == g2.genes.end())
				disjoint_count++;

		for (; it2 != g2.genes.end(); it2++)
			if (g1.genes.find((*it2).second.innovation_num) == g1.genes.end())
				disjoint_count++;

		return (1. * disjoint_count) /
			   (1. * std::max(g1.genes.size(), g2.genes.size()));
	}

	double pool::weights(const genome &g1, const genome &g2)
	{
		auto it1 = g1.genes.begin();

		double sum = 0.0;
		unsigned int coincident = 0;

		for (; it1 != g1.genes.end(); it1++)
		{
			auto it2 = g2.genes.find((*it1).second.innovation_num);
			if (it2 != g2.genes.end())
			{
				coincident++;
				sum += std::abs((*it1).second.weight - (*it2).second.weight);
			}
		}

		return 1. * sum / (1. * coincident);
	}

	bool pool::is_same_species(const genome &g1, const genome &g2)
	{
		double dd = this->speciating_parameters.delta_disjoint * disjoint(g1, g2);
		double dw = this->speciating_parameters.delta_weights * weights(g1, g2);
		return dd + dw < this->speciating_parameters.delta_threshold;
	}

	void pool::rank_globally()
	{
		std::vector<genome *> global;
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			for (size_t i = 0; i < (*s).genomes.size(); i++)
				global.push_back(&((*s).genomes[i]));

		std::sort(global.begin(), global.end(), [](genome *&a, genome *&b) -> bool
				  { return a->fitness < b->fitness; });
		for (size_t j = 0; j < global.size(); j++)
			global[j]->global_rank = j + 1;
	}

	void pool::calculate_average_fitness(specie &s)
	{
		unsigned int total = 0;
		for (size_t i = 0; i < s.genomes.size(); i++)
			total += s.genomes[i].global_rank;
		s.average_fitness = total / s.genomes.size();
	}

	unsigned int pool::total_average_fitness()
	{
		unsigned int total = 0;
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			total += (*s).average_fitness;
		return total;
	}

	void pool::cull_species(bool cut_to_one)
	{
		for (auto s = this->species.begin(); s != this->species.end(); s++)
		{
			std::sort((*s).genomes.begin(), (*s).genomes.end(),
					  [](genome &a, genome &b)
					  { return a.fitness > b.fitness; });

			unsigned int remaining = std::ceil((*s).genomes.size() * 1.0 / 2.0);
			// this will leave the most fit genome in specie,
			// letting him make more and more babies (until someone in
			// specie beat him or he becomes weaker during mutations
			if (cut_to_one)
				remaining = 1;
			while ((*s).genomes.size() > remaining)
				(*s).genomes.pop_back();
		}
	}

	genome pool::breed_child(specie &s)
	{
		genome child(this->network_info, this->mutation_rates);
		std::uniform_real_distribution<double> distributor(0.0, 1.0);
		std::uniform_int_distribution<unsigned int> choose_genome(
			0, s.genomes.size() - 1);
		if (distributor(this->generator) < this->mutation_rates.crossover_chance)
		{
			unsigned int g1id, g2id;
			genome &g1 = s.genomes[g1id = choose_genome(this->generator)];
			genome &g2 = s.genomes[g2id = choose_genome(this->generator)];

			// QUESTION: if g1 == g2, then you can make a baby by fapping?
			child = this->crossover(g1, g2);
		}
		else
		{
			genome &g = s.genomes[choose_genome(this->generator)];
			child = g;
		}

		this->mutate(child);
		return child;
	}

	void pool::remove_stale_species()
	{
		auto s = this->species.begin();
		while (s != this->species.end())
		{
			genome &g = *(std::max_element(
				(*s).genomes.begin(), (*s).genomes.end(),
				[](genome &a, genome &b) -> bool
				{ return a.fitness < b.fitness; }));

			if (g.fitness > (*s).top_fitness)
			{
				(*s).top_fitness = g.fitness;
				(*s).staleness = 0;
			}
			else
				(*s).staleness++;

			if (!((*s).staleness < this->speciating_parameters.stale_species ||
				  (*s).top_fitness >= this->max_fitness))
				this->species.erase(s++);
			else
				s++;
		}
	}

	void pool::remove_weak_species()
	{
		unsigned int sum = this->total_average_fitness();
		auto s = this->species.begin();
		while (s != this->species.end())
		{
			double breed = std::floor((1. * (*s).average_fitness) / (1. * sum) * 1. *
									  this->speciating_parameters.population);
			if (breed >= 1.0)
				s++;
			else
				this->species.erase(s++);
		}
	}

	void pool::add_to_species(genome &child)
	{
		auto s = this->species.begin();
		while (s != this->species.end())
		{
			if (this->is_same_species(child, (*s).genomes[0]))
			{
				(*s).genomes.push_back(child);
				break;
			}
			++s;
		}

		if (s == this->species.end())
		{
			specie new_specie;
			new_specie.genomes.push_back(child);
			this->species.push_back(new_specie);
		}
	}

	void pool::new_generation()
	{

		this->innovation.reset();
		this->cull_species(false);
		this->rank_globally();
		this->remove_stale_species();

		for (auto s = this->species.begin(); s != this->species.end(); s++)
			this->calculate_average_fitness((*s));
		this->remove_weak_species();

		std::vector<genome> children;
		unsigned int sum = this->total_average_fitness();
		for (auto s = this->species.begin(); s != this->species.end(); s++)
		{
			unsigned int breed =
				std::floor(((1. * (*s).average_fitness) / (1. * sum)) * 1. *
						   this->speciating_parameters.population) -
				1;
			for (unsigned int i = 0; i < breed; i++)
				children.push_back(this->breed_child(*s));
		}

		this->cull_species(true); // now in each species we have only one genome

		// preparing for MAKING BABIES <3
		std::uniform_int_distribution<unsigned int> choose_specie(
			0, this->species.size() - 1);
		std::vector<specie *> species_pointer(0);
		for (auto s = this->species.begin(); s != this->species.end(); s++)
			species_pointer.push_back(&(*s));
		if (this->species.size() == 0)
			printf("[ERROR] Wtf? Zero species in the world! All dead? Where is that "
				   "fucking NOAH and his fucking boat?\n");
		else
			while (children.size() + this->species.size() <
				   this->speciating_parameters.population)
				children.push_back(
					this->breed_child(*species_pointer[choose_specie(this->generator)]));

		for (size_t i = 0; i < children.size(); i++)
			this->add_to_species(children[i]);
		this->generation_number++;
	}
#if INDEPENDENT
	void pool::import_fromfile(std::virtual_file &file)
	{
		file.rewind();
		this->species.clear();
		// current state
		unsigned int innovation_num;
		file.scanf("%u", &innovation_num);
		this->innovation.set_innovation_number(innovation_num);
		file.scanf("%u%u", &this->generation_number, &this->max_fitness);

		// network information
		file.scanf("%u%u%u", &this->network_info.input_size,
				   &this->network_info.output_size, &this->network_info.bias_size);
		this->network_info.functional_nodes = this->network_info.input_size +
											  this->network_info.output_size +
											  this->network_info.bias_size;

		std::string rec = file.get_str();
		if (rec == "rec")
			this->network_info.recurrent = true;
		if (rec == "nonrec")
			this->network_info.recurrent = false;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		if (file.get_str() != "CHANGEABLE_ACTIVATION_AND_AGGREGATION")
		{
			printf("[ERROR] Failed to import a neat::pool without CHANGEABLE_ACTIVATION_AND_AGGREGATION enabled");
			exit(-1);
		}
#endif
		// population information
		this->speciating_parameters.read(file);

		// mutation parameters
		this->mutation_rates.read(file);

		// species information
		unsigned int species_number;
		file.scanf("%u", &species_number);
		this->species.clear();

		for (unsigned int c = 0; c < species_number; c++)
		{
			specie new_specie;
#ifdef GIVING_NAMES_FOR_SPECIES
			new_specie.name = file.get_str();
#endif

			unsigned int specie_population;

			file.scanf("%u%u%u%u", &new_specie.top_fitness, &new_specie.average_fitness,
					   &new_specie.staleness, &specie_population);

			for (unsigned int i = 0; i < specie_population; i++)
			{
				genome new_genome(this->network_info, this->mutation_rates);

				file.scanf("%u%u%u", &new_genome.fitness, &new_genome.adjusted_fitness,
						   &new_genome.global_rank);

				new_genome.mutation_rates.read(file);

				unsigned int gene_number;
				file.scanf("%u", &gene_number);

				for (unsigned int j = 0; j < gene_number; j++)
				{
					gene new_gene;

					file.scanf("%u%u%u%lf%c", &new_gene.innovation_num, &new_gene.from_node,
							   &new_gene.to_node, &new_gene.weight, &new_gene.enabled);
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
					new_gene.aggregation_name = file.get_str();
					new_gene.activation_name = file.get_str();
					new_gene.aggregation = this->aggregation_funcs.at(new_gene.aggregation_name);
					new_gene.activation = this->activation_funcs.at(new_gene.activation_name);
#endif
					new_genome.genes[new_gene.innovation_num] = new_gene;
				}

				new_specie.genomes.push_back(new_genome);
			}

			this->species.push_back(new_specie);
		}
	}
#else
	void pool::import_fromfile(std::string filename)
	{
		std::ifstream input;
		input.open(filename);
		if (!input.is_open())
		{
			std::cerr << "cannot open file '" << filename << "' !\n";
			return;
		}

		this->species.clear();
		try
		{
			// current state
			unsigned int innovation_num;
			input >> innovation_num;
			this->innovation.set_innovation_number(innovation_num);
			input >> this->generation_number;
			input >> this->max_fitness;

			// network information
			input >> this->network_info.input_size >> this->network_info.output_size >>
				this->network_info.bias_size;
			this->network_info.functional_nodes = this->network_info.input_size +
												  this->network_info.output_size +
												  this->network_info.bias_size;

			std::string rec;
			input >> rec;
			if (rec == "rec")
				this->network_info.recurrent = true;
			if (rec == "nonrec")
				this->network_info.recurrent = false;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
			input >> rec;
			if (rec != "CHANGEABLE_ACTIVATION_AND_AGGREGATION")
			{
				throw "[ERROR] Failed to import a neat::pool without CHANGEABLE_ACTIVATION_AND_AGGREGATION enabled\n";
			}
#endif
			// population information
			this->speciating_parameters.read(input);

			// mutation parameters
			this->mutation_rates.read(input);

			// species information
			unsigned int species_number;
			input >> species_number;
			this->species.clear();

			for (unsigned int c = 0; c < species_number; c++)
			{
				specie new_specie;
#ifdef GIVING_NAMES_FOR_SPECIES
				input >> new_specie.name;
#endif
				input >> new_specie.top_fitness;
				input >> new_specie.average_fitness;
				input >> new_specie.staleness;

				unsigned int specie_population;
				input >> specie_population;

				for (unsigned int i = 0; i < specie_population; i++)
				{
					genome new_genome(this->network_info, this->mutation_rates);
					input >> new_genome.fitness;
					input >> new_genome.adjusted_fitness;
					input >> new_genome.global_rank;

					new_genome.mutation_rates.read(input);

					unsigned int gene_number;
					input >> new_genome.max_neuron >> gene_number;

					for (unsigned int j = 0; j < gene_number; j++)
					{
						gene new_gene;
						input >> new_gene.innovation_num;
						input >> new_gene.from_node;
						input >> new_gene.to_node;
						input >> new_gene.weight;
						input >> new_gene.enabled;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
						input >> new_gene.aggregation_name;
						input >> new_gene.activation_name;
						new_gene.aggregation = this->aggregation_funcs.at(new_gene.aggregation_name);
						new_gene.activation = this->activation_funcs.at(new_gene.activation_name);
#endif
						new_genome.genes[new_gene.innovation_num] = new_gene;
					}

					new_specie.genomes.push_back(new_genome);
				}

				this->species.push_back(new_specie);
			}
		}
		catch (std::string error_message)
		{
			std::cerr << error_message;
		}

		input.close();
	}
#endif
#if INDEPENDENT
	void pool::export_tofile(std::virtual_file &file)
	{

		// current state
		file.concat(std::to_string(this->innovation.number()) + "\n" +
					std::to_string(this->generation_number) + "\n" +
					std::to_string(this->max_fitness) + "\n");

		// network information
		file.concat(std::to_string(this->network_info.input_size) + " " +
					std::to_string(this->network_info.output_size) + " " +
					std::to_string(this->network_info.bias_size) + "\n");
		if (this->network_info.recurrent)
			file.concat("rec\n");
		else
			file.concat("nonrec\n");
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		file.concat("CHANGEABLE_ACTIVATION_AND_AGGREGATION\n");
#endif
		this->network_info.functional_nodes = this->network_info.input_size +
											  this->network_info.output_size +
											  this->network_info.bias_size;

		// population information
		this->speciating_parameters.write(file, "");

		// mutation parameters
		this->mutation_rates.write(file, "");

		// species information
		file.concat(std::to_string(this->species.size()) + "\n");
		for (auto s = this->species.begin(); s != this->species.end(); s++)
		{
			file.concat("   ");
#ifdef GIVING_NAMES_FOR_SPECIES
			// why not use the arrow?
			// leaving this here anyway

			// bugfix: this needs to end with \n not a space for it to be read correctly
			file.concat((*s).name + "\n");
#endif
			file.concat(std::to_string((*s).top_fitness) + " " +
						std::to_string((*s).average_fitness) + " " +
						std::to_string((*s).staleness) + "\n");

			file.concat("   " + std::to_string((*s).genomes.size()) + "\n");
			for (size_t i = 0; i < (*s).genomes.size(); i++)
			{
				file.concat("      " + std::to_string((*s).genomes[i].fitness) + " " +
							std::to_string((*s).genomes[i].adjusted_fitness) + " " +
							std::to_string((*s).genomes[i].global_rank) + "\n");

				(*s).genomes[i].mutation_rates.write(file, "      ");

				file.concat("      " + std::to_string((*s).genomes[i].max_neuron) + " " +
							std::to_string((*s).genomes[i].genes.size()) + "\n");
				for (auto it = (*s).genomes[i].genes.begin();
					 it != (*s).genomes[i].genes.end(); it++)
				{
					gene &g = (*it).second;
					file.concat("         " + std::to_string(g.innovation_num) + " " +
								std::to_string(g.from_node) + " " +
								std::to_string(g.to_node) + " " + std::to_string(g.weight) +
								" " + std::to_string(g.enabled) + "\n");
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
					file.concat(g.aggregation_name + "\n" + g.activation_name + "\n");
#endif
				}
			}

			file.concat("\n\n");
		}
	}
#else
	void pool::export_tofile(std::string filename)
	{
		std::ofstream output;
		output.open(filename);
		if (!output.is_open())
		{
			std::cerr << "cannot open file '" << filename << "' !";
			return;
		}

		// current state
		output << this->innovation.number() << "\n";
		output << this->generation_number << "\n";
		output << this->max_fitness << "\n";

		// network information
		output << this->network_info.input_size << " "
			   << this->network_info.output_size << " "
			   << this->network_info.bias_size << "\n";
		if (this->network_info.recurrent)
			output << "rec"
				   << "\n";
		else
			output << "nonrec"
				   << "\n";
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		output << "CHANGEABLE_ACTIVATION_AND_AGGREGATION\n";
#endif
		this->network_info.functional_nodes = this->network_info.input_size +
											  this->network_info.output_size +
											  this->network_info.bias_size;

		// population information
		this->speciating_parameters.write(output, "");

		// mutation parameters
		this->mutation_rates.write(output, "");

		// species information
		output << this->species.size() << "\n";
		for (auto s = this->species.begin(); s != this->species.end(); s++)
		{
			output << "   ";
#ifdef GIVING_NAMES_FOR_SPECIES
			output << (*s).name << "\n";
#endif
			output << (*s).top_fitness << " ";
			output << (*s).average_fitness << " ";
			output << (*s).staleness << "\n";

			output << "   " << (*s).genomes.size() << "\n";
			for (size_t i = 0; i < (*s).genomes.size(); i++)
			{
				output << "      ";
				output << (*s).genomes[i].fitness << " ";
				output << (*s).genomes[i].adjusted_fitness << " ";
				output << (*s).genomes[i].global_rank << "\n";

				(*s).genomes[i].mutation_rates.write(output, "      ");

				output << "      " << (*s).genomes[i].max_neuron << " "
					   << (*s).genomes[i].genes.size() << "\n";
				for (auto it = (*s).genomes[i].genes.begin();
					 it != (*s).genomes[i].genes.end(); it++)
				{
					gene &g = (*it).second;
					output << "         ";
					output << g.innovation_num << " " << g.from_node << " " << g.to_node
						   << " " << g.weight << " " << g.enabled << "\n";
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
					output << g.aggregation_name << "\n"
						   << g.activation_name << "\n";
#endif
				}
			}

			output << "\n"
				   << "\n";
		}
		output.close();
	}
#endif
// DO NOT rewind for these
#if INDEPENDENT
	void mutation_rate_container::read(std::virtual_file &file)
	{
		file.scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf", &this->connection_mutate_chance,
				   &this->perturb_chance, &this->crossover_chance,
				   &this->link_mutation_chance, &this->node_mutation_chance,
				   &this->bias_mutation_chance, &this->step_size,
				   &this->disable_mutation_chance, &this->enable_mutation_chance);
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		file.scanf("%lf%lf", &this->activation_mutate_chance, &this->aggregation_mutate_chance);
#endif
	}
#else
	void mutation_rate_container::read(std::ifstream &o)
	{
		o >> this->connection_mutate_chance;
		o >> this->perturb_chance;
		o >> this->crossover_chance;
		o >> this->link_mutation_chance;
		o >> this->node_mutation_chance;
		o >> this->bias_mutation_chance;
		o >> this->step_size;
		o >> this->disable_mutation_chance;
		o >> this->enable_mutation_chance;
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		o >> this->activation_mutate_chance;
		o >> this->aggregation_mutate_chance;
#endif
	}
#endif

#if INDEPENDENT
	void mutation_rate_container::write(std::virtual_file &file,
										std::string prefix)
	{
		file.concat(prefix + std::to_string(this->connection_mutate_chance) + "\n" +
					prefix + std::to_string(this->perturb_chance) + "\n" + prefix +
					std::to_string(this->crossover_chance) + "\n" + prefix +
					std::to_string(this->link_mutation_chance) + "\n" + prefix +
					std::to_string(this->node_mutation_chance) + "\n" + prefix +
					std::to_string(this->bias_mutation_chance) + "\n" + prefix +
					std::to_string(this->step_size) + "\n" + prefix +
					std::to_string(this->disable_mutation_chance) + "\n" + prefix +
					std::to_string(this->enable_mutation_chance) + "\n");
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		file.concat(prefix + std::to_string(this->activation_mutate_chance) + "\n" +
					prefix + std::to_string(this->aggregation_mutate_chance) + "\n");
#endif
	}
#else
	void mutation_rate_container::write(std::ofstream &o, std::string prefix)
	{
		o << prefix << this->connection_mutate_chance << "\n";
		o << prefix << this->perturb_chance << "\n";
		o << prefix << this->crossover_chance << "\n";
		o << prefix << this->link_mutation_chance << "\n";
		o << prefix << this->node_mutation_chance << "\n";
		o << prefix << this->bias_mutation_chance << "\n";
		o << prefix << this->step_size << "\n";
		o << prefix << this->disable_mutation_chance << "\n";
		o << prefix << this->enable_mutation_chance << "\n";
#ifdef CHANGEABLE_ACTIVATION_AND_AGGREGATION
		o << prefix << this->activation_mutate_chance << "\n";
		o << prefix << this->aggregation_mutate_chance << "\n";
#endif
	}
#endif

#if INDEPENDENT
	void speciating_parameter_container::read(std::virtual_file &file)
	{
		file.scanf("%u%lf%lf%lf%u", &this->population, &this->delta_disjoint,
				   &this->delta_weights, &this->delta_threshold,
				   &this->stale_species);
	}
#else
	void speciating_parameter_container::read(std::ifstream &o)
	{
		o >> this->population;
		o >> this->delta_disjoint;
		o >> this->delta_weights;
		o >> this->delta_threshold;
		o >> this->stale_species;
	}
#endif

#if INDEPENDENT
	void speciating_parameter_container::write(std::virtual_file &file,
											   std::string prefix)
	{
		file.concat(prefix + std::to_string(this->population) + "\n" + prefix +
					std::to_string(this->delta_disjoint) + "\n" + prefix +
					std::to_string(this->delta_weights) + "\n" + prefix +
					std::to_string(this->delta_threshold) + "\n" + prefix +
					std::to_string(this->stale_species) + "\n");
	}
#else
	void speciating_parameter_container::write(std::ofstream &o,
											   std::string prefix)
	{
		o << prefix << this->population << "\n";
		o << prefix << this->delta_disjoint << "\n";
		o << prefix << this->delta_weights << "\n";
		o << prefix << this->delta_threshold << "\n";
		o << prefix << this->stale_species << "\n";
	}
#endif

} // end of namespace neat

#endif
