/*
 * Node.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */

#include "Node.h"

Node::Node() {
	// TODO Auto-generated constructor stub.
	this->out_degree = 0;
}

Node::~Node() {
	// TODO Auto-generated destructor stub.
}

void Node::addIncomingEdge(int incoming) {
	this->in_nodes.push_back(incoming);
}

std::vector<int> Node::getIncomingNodes() {
	return this->in_nodes;
}

void Node::incrementOutDegree() {
	++(this->out_degree);
}

int Node::getOutDegree() {
	return out_degree;
}




