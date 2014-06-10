/*
 * AdjacencyListGraph.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: gnanda
 */
#include <cstdlib>

#include "AdjacencyListGraph.h"

AdjacencyListGraph::AdjacencyListGraph() {
	// TODO Auto-generated constructor stub
	count = 0;
}

AdjacencyListGraph::~AdjacencyListGraph() {
	// TODO Auto-generated destructor stub
}

// Read the graph input file and generate a map from node number
// to the actual node objects.
void AdjacencyListGraph::ingestGraph(std::string path) {
	ifstream myfile (path.c_str());

	string line;
	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			string delimiter = " ";

			// Get start and end of an edge.
			int start = atol(line.substr(0, line.find(delimiter)).c_str());
			int end = atol(line.substr(line.find(delimiter) + 1, line.size()).c_str());

			Node* edgeStart = NULL;
			Node* edgeEnd   = NULL;

			// If the starting node does not exist, create it
			// and add it to vector.
			if (start + 1 > nodeVector.size()) {
				nodeVector.resize(start + 1, NULL);
			}

			if (nodeVector[start] == NULL) {
				edgeStart = new Node();
				this->nodeVector[start] = edgeStart;
				++count;
			}

			if (end + 1 > nodeVector.size()) {
				nodeVector.resize(end + 1, NULL);
			}

			if (nodeVector[end] == NULL) {
				edgeEnd = new Node();
				this->nodeVector[end] = edgeEnd;
				++count;
			}

//			if (!nodeMap.count(start)) {
//				edgeStart = new Node();
//				this->nodeMap[start] = edgeStart;
//			}
//
//			// Now create the ending node and add it to the map.
//			if (!nodeMap.count(end)) {
//				edgeEnd = new Node();
//				this->nodeMap[end] = edgeEnd;
//			}

//			edgeStart = this->nodeMap[start];
//			edgeEnd = this->nodeMap[end];

			edgeStart = this->nodeVector[start];
			edgeEnd = this->nodeVector[end];

			edgeStart->incrementOutDegree();
			edgeEnd->addIncomingEdge(start);
		}
		myfile.close();
	}

	else {
		cout << "Unable to open: " << path << endl;
	}
}

Node* AdjacencyListGraph::getNode(int number) {
	return this->nodeVector[number];
}

int AdjacencyListGraph::getTotalNodes() {
	return count;
}

int AdjacencyListGraph::getVectorLength() {
	return this->nodeVector.size();
}
