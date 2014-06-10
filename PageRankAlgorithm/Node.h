/*
 * Node.h
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */
#include <vector>

#ifndef NODE_H_
#define NODE_H_

/**
 * This class represents a node of the graph.
 */
class Node {
private:
	std::vector<int> in_nodes; // Set of nodes which have incoming nodes to this node.
	int out_degree; // Represents the count of outgoing nodes.
public:
	Node();
	virtual ~Node();

	// Adds an incoming edge to this node.
	void addIncomingEdge(int incoming);

	// Returns a vector of the incoming nodes.
	std::vector<int> getIncomingNodes();

	// Increment the out degree.
	void incrementOutDegree();

	// Returns the outDegree.
	int getOutDegree();
};

#endif /* NODE_H_ */
