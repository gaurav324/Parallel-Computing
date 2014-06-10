/*
 * AdjacencyListGraph.h
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */

#ifndef ADJACENCYLISTGRAPH_H_
#define ADJACENCYLISTGRAPH_H_

#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "Graph.h"
#include "Node.h"

using namespace std;

class AdjacencyListGraph: public Graph {
private:
	// Vector of all the nodes present in the graph.
	vector<Node*> nodeVector;

	// Total nodes.
	int count;
public:
	AdjacencyListGraph();
	virtual ~AdjacencyListGraph();

	// Read a graph and load in memory.
	void ingestGraph(string path);

	// Given the node index, it returns the Node*.
	Node* getNode(int index);

	// Returns the number of nodes.
	int getTotalNodes();

	// Return the length of vector which might also be having some empty spaces.
	int getVectorLength();
};

#endif /* ADJACENCYLISTGRAPH_H_ */
