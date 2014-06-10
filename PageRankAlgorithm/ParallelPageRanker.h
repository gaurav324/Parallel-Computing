/*
 * ParallelPageRanker.h
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */
#include <iostream>
#include <vector>
#include <algorithm>

#include <omp.h>
#include "tgmath.h"
#include "time.h"
#include "assert.h"
#include "AdjacencyListGraph.h"

#ifndef PARALLELPAGERANKER_H_
#define PARLLELPAGERANKER_H_

#define CONVERGENCE 1.0 / 10000
#define DAMPING_FACTOR 0.85
#define MAX_ITERATIONS 100

using namespace std;

class ParallelPageRanker {
private:
	AdjacencyListGraph graph;

	// Store ranks for all the nodes.
	vector<float>* rank;

	// Reads the graph and populates the rank matrix.
	void generateRankMatrix();
public:
	ParallelPageRanker();
	virtual ~ParallelPageRanker();

	// Attach a graph file to the PageRanker.
	void attachGraphFile(string File);

	// This function would time the pageRank Convergence.
	map<int, float> timePageRankTwoCopy();

	// This function would time the pageRank Convergence.
	void timePageRankOneCopy(map<int, float>& twoCopyTop);
};

#endif /* SERIALPAGERANKER_H_ */
