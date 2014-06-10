#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <queue>

#include "omp.h"
#include "limits.h"

using namespace std;

class end
{
public:
	int dest;
	int weight;

	end(int dest, int weight) {
		this->dest = dest;
		this->weight = weight;
	}
};

typedef pair<int,int> wn;

class graph
{
public:
	// Total vertices and edges.
	int total_vertices;
	int total_edges;

	vector< vector<end*>* > vertices;

	graph(int vertices, int edge_count) {
		total_vertices = vertices;
		total_edges = edge_count;
	}

	void add(int source, int dest, int weight) {
		if (source >= vertices.size()) {
			vertices.resize(source + 1, NULL);
		}

		if (dest >= vertices.size()) {
			vertices.resize(dest + 1, NULL);
		}

		end* e = new end(dest, weight);

		if (vertices[source] == NULL) {
			vertices[source] = new vector<end*>();
		}
		vertices[source]->push_back(e);
	}
};

int main(int argc, char* argv[]) {

	int source_node = atol(argv[2]);
	string filename = argv[3];

	cout << "Source node: " << source_node << endl;
	cout << "filename: " << filename << endl;

	ifstream myfile (argv[1]);

	string line;
	string delimiter = " ";

	cout << "Reading in graph." << endl;
	if (myfile.is_open()) {
		getline(myfile, line);
		int vertices = atol(line.substr(0, line.find(delimiter)).c_str());
		int edge_count = atol(line.substr(line.find(delimiter) + 1, line.size()).c_str());

		graph G(vertices, edge_count);

		while (getline(myfile, line)) {
			// Get start, end and weight of an edge.
			int source = atol(line.substr(0, line.find(delimiter)).c_str());
			line = line.substr(line.find(delimiter) + 1, line.size());

			int end = atol(line.substr(0, line.find(delimiter)).c_str());

			int weight = atol(line.substr(line.find(delimiter) + 1, line.size()).c_str());

			G.add(source, end, weight);
		}
		myfile.close();

		cout << "Total vertices: " << G.total_vertices << endl;
		cout << "Total edges: " << G.total_edges << endl;

		int dist[G.total_vertices];
		for (int i=0; i < G.total_vertices;++i) {
			dist[i] = INT_MAX;
		}

		double start_time;
		start_time = omp_get_wtime();

		cout << "Dijkstra starts here. " << endl;

		priority_queue<wn, std::vector<wn>, std::greater<wn> > Q;
		Q.push(wn(0, source_node));

		dist[source_node] = 0;
		while(!Q.empty()) {
			wn weight_node = Q.top();
			Q.pop();

			int weight = weight_node.first;
			int src = weight_node.second;

			if(weight <= dist[src]) {
				vector<end*>* ends = G.vertices[src];
				if (ends == NULL || ends->size() == 0){
					continue;
				}
				for (int i=0; i < ends->size(); ++i) {
					end* e = (*ends)[i];

					int new_weight = dist[src] + e->weight;
					if (new_weight < dist[e->dest]) {
						dist[e->dest] = new_weight;
						Q.push(wn(dist[e->dest], e->dest));
					}
				}
			}
		}

		double end_time;
		end_time = omp_get_wtime();
		// BELLMAN FORD ENDS HERE.

		cout << "Time taken for serial(sec): " << end_time - start_time << endl;

		ofstream myfile;
		myfile.open (filename.c_str());
		for (int i=0; i < G.total_vertices; i++) {
			myfile << i << "," << dist[i] << endl;
		}
		myfile.close();
	}

	return 1;
}
