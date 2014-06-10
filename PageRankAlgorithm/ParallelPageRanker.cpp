/*
 * ParallelPageRanker.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */

#include "ParallelPageRanker.h"

ParallelPageRanker::ParallelPageRanker() {
	// TODO Auto-generated constructor stub
	std::cout << "Constructing Parallel Page Ranker." << endl;
	rank = new vector<float>();
}

ParallelPageRanker::~ParallelPageRanker() {
	// TODO Auto-generated destructor stub
	std::cout << "Destructing Parallel Page Ranker." << endl;
}

void ParallelPageRanker::attachGraphFile(string file) {
	//std::cout << "ParallelPageRanker: Going to ingestGraph." << endl;
	this->graph.ingestGraph(file);
	//std::cout << "ParallelPageRanker: Done with ingestGraph." << endl;

	// Ensure that the graph read the file properly.
	//std::cout << "ParallelPageRanker: Going to assert." << endl;
	assert(this->graph.getTotalNodes() > 0);
	//std::cout << "ParallelPageRanker: Assertion Complete." << endl;

	// Let us also generate the rank matrix here itself.
	//std::cout << "ParallelPageRanker: Going to generateRankMatrix." << endl;
	generateRankMatrix();
	//std::cout << "ParallelPageRanker: generateRankMatrix Complete." << endl;
}

void ParallelPageRanker::generateRankMatrix() {
	float initial = (1.0 / this->graph.getTotalNodes());
	for (int l=0; l < this->graph.getVectorLength(); ++l) {
		this->rank->push_back(initial);
	}
}

map<int, float> ParallelPageRanker::timePageRankTwoCopy() {
	double start_time;
	start_time = omp_get_wtime();

	float initial = 1.0 / this->graph.getTotalNodes();

	// Store ranks(i+1) here.
	vector<float>* rank1 = new vector<float>();
	cout << "Total Nodes are: " << this->graph.getTotalNodes() << endl;
	for (int l=0; l < this->graph.getVectorLength(); l++) {
			(*this->rank)[l] = initial;
			rank1->push_back(0.0);
	}

	// Compute the constant factor and store in advance,
	float const_factor = (1.0 - DAMPING_FACTOR) / this->graph.getTotalNodes();

	vector<float> diffs;
	diffs.resize(this->graph.getVectorLength());

	for(int i=0; i < MAX_ITERATIONS; ++i) {
		float total_change = 0.0;

		// Iterate over all the members in the rank vector.
		#pragma omp parallel for num_threads(60)
		for (int l=0; l < (int)this->rank->size(); ++l) {
			// Get the list of the neighbors in context.
			Node* x = this->graph.getNode(l);
			if (x == NULL) {
				(*rank)[l] = 0.0;
				//cout << "No node present for number: " << l;
				continue;
			}

			std::vector<int> in_nodes = x->getIncomingNodes();

			// Sum up the rank-factor of all the incoming nodes.
			float moving_factor = 0.0;
			for(std::vector<int>::iterator it=in_nodes.begin();
					it != in_nodes.end(); ++it) {
				Node* incoming_node = graph.getNode(*it);
				moving_factor += ((*rank)[*it] / incoming_node->getOutDegree());
			}

			float new_rank = const_factor +
					DAMPING_FACTOR * moving_factor;
			(*rank1)[l] = new_rank;
			float change = new_rank - (*rank)[l];
			diffs[l] = change;
		}

		for (int l=0; l < (int)this->rank->size(); ++l) {
			total_change += diffs[l] > 0 ? diffs[l] : diffs[l] * -1;
		}

		cout << "Total change: " << total_change << endl;
 		if (total_change < (CONVERGENCE)) {
			break;
		}

		vector<float>* temp = rank1;
		rank1 = rank;
		rank = temp;
	}

	cout << "Time taken: "
			<< (float) (omp_get_wtime() - start_time)
			<< " sec" << endl;

	// This would sort the vector. But would also retain the indices.
	std::map<float, int> twoCopy;
	for (int l=0; l < this->graph.getVectorLength(); ++l) {
		twoCopy[(*rank1)[l]] = l;
	}

	// Return the top 100 page ranks.
	std::map<int, float> result;
	int counter = 0;
	for (std::map<float, int>::reverse_iterator rit=twoCopy.rbegin();
			rit!=twoCopy.rend(); ++rit) {
		result[rit->second] = rit->first;
		//cout << rit->first << endl;
		++counter;
		if (counter == 100) {
			break;
		}
	}
	return result;
}

void ParallelPageRanker::timePageRankOneCopy(map<int, float>& twoCopy) {
	float initial = 1.0 / this->graph.getTotalNodes();
	for (int l=0; l < this->graph.getVectorLength(); l++) {
		(*this->rank)[l] = initial;
	}

	double start_time;
	start_time = omp_get_wtime();

	// Compute the constant factor and store in advance,
	float const_factor = (1 - DAMPING_FACTOR) / this->graph.getTotalNodes();
	cout << "Total Nodes are: " << this->graph.getTotalNodes() << endl;

	// Iterate over all the members in the rank vector.
	vector<float> diffs;
	diffs.resize(this->graph.getVectorLength());

	for(int i=0; i < MAX_ITERATIONS; ++i) {
		float total_change = 0.0;

		#pragma omp parallel for num_threads(60)
		for (int l=0; l < (int)this->rank->size(); ++l) {
			// Get the list of the neighbors in context.
			Node* x = this->graph.getNode(l);
			if (x == NULL) {
				(*rank)[l] = 0.0;
				continue;
			}

			std::vector<int> in_nodes = x->getIncomingNodes();

			// Sum up the rank-factor of all the incoming nodes.
			float moving_factor = 0.0;
			for(std::vector<int>::iterator it=in_nodes.begin();
					it != in_nodes.end(); ++it) {
				Node* incoming_node = graph.getNode(*it);
				moving_factor += ((*rank)[*it] / incoming_node->getOutDegree());
			}

			float new_rank = const_factor +
					DAMPING_FACTOR * moving_factor;

			float change = new_rank - (*rank)[l];
			diffs[l] = change;

			(*rank)[l] = new_rank;
		}

		for (int l=0; l < (int)this->rank->size(); ++l) {
			total_change += diffs[l] > 0 ? diffs[l] : diffs[l] * -1;
		}

		cout << "Total change: " << total_change << endl;
 		if (total_change < (CONVERGENCE)) {
			break;
		}
	}


	cout << "Time taken: "
		 << (float) (omp_get_wtime() - start_time)
		 << " sec" << endl;

	// Find Kendaull's Tau Distance.
	int count = 0;
	float distance = 0.0;
	for (map<int, float>::iterator it1=twoCopy.begin(); it1!=twoCopy.end(); ++it1) {
		for (map<int, float>::iterator it2=twoCopy.begin(); it2!=it1; ++it2) {
			if ((it1->second < it2->second) && ((*rank)[it1->first] > (*rank)[it2->first])) {
				//cout << "Two copy values: "<< it1->second << " " << it2->second << endl;
				//cout << "One copy values: " << (*rank)[it1->first] << " " << (*rank)[it2->first] << endl;
				++count;
			}
			else if ((it1->second > it2->second) && ((*rank)[it1->first] < (*rank)[it2->first])) {
				//cout << "Two copy values: "<< it1->second << " " << it2->second << endl;
				//cout << "One copy values: " << (*rank)[it1->first] << " " << (*rank)[it2->first] << endl;
				++count;
			}
		}
	}

	distance = count * 1.0 / (50 * 99);
	cout << "Kendaull's Tau Distance: " << distance << endl;
}
