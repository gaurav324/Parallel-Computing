#include <iostream>

// Local Imports.
#include "SerialPageRanker.h"
#include "ParallelPageRanker.h"

#include "AdjacencyListGraph.h"
int main(int argc, char* argv[]) {

	if (argc > 2) {
		ParallelPageRanker ppr;

		cout << "Attaching graph file." << endl;
		ppr.attachGraphFile(argv[1]);

		cout << "Running Two page Rank with two copies of Rank Matrix:" << endl;
		// Returns the node node number against, page rank hash_map.
		// This hash_map has only top 100 elements.
		map<int, float> twoCopy = ppr.timePageRankTwoCopy();

		cout << "Running Two page Rank with one copy of Rank Matrix:" << endl;
		ppr.timePageRankOneCopy(twoCopy);
	}

	else {
		SerialPageRanker spr;

		cout << "Attaching graph file." << endl;
		spr.attachGraphFile(argv[1]);

		cout << "Running Two page Rank with two copies of Rank Matrix:" << endl;
		// Returns the node node number against, page rank hash_map.
		// This hash_map has only top 100 elements.
		map<int, float> twoCopy = spr.timePageRankTwoCopy();

		cout << "Running Two page Rank with one copy of Rank Matrix:" << endl;
		spr.timePageRankOneCopy(twoCopy);
	}
	return 0;
}
