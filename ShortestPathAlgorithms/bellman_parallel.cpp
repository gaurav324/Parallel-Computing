#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>
#include <ctime>

#include "omp.h"
#include "limits.h"
using namespace std;


class edge
{
public:
	int src;
	int dest;
	int weight;

	edge(int source, int vertex, int weight) {
		src = source;
		dest = vertex;
		this->weight = weight;
	}
};

class graph
{
public:
	// Total vertices and edges.
	int total_vertices;
	int total_edges;

	// Array of the edges.
	vector<edge*> edges;

	graph(int vertices, int edge_count) {
		total_vertices = vertices;
		total_edges = edge_count;
	}

	void add(int source, int vertex, int weight) {
		edge *ed = new edge(source, vertex, weight);
		edges.push_back(ed);
	}
};
int main(int argc, char* argv[]) {

    int source_node = atol(argv[3]);
	string filename = argv[4];

    cout << "Source node: " << source_node << endl;
	cout << "filename: " << filename << endl;

	ifstream myfile (argv[1]);

	string line;
	string delimiter = " ";

	if (myfile.is_open()) {
		getline(myfile, line);
		int vertices = atol(line.substr(0, line.find(delimiter)).c_str());
		int edge_count = atol(line.substr(line.find(delimiter) + 1, line.size()).c_str());

		// READ THE GRAPH.
        cout << "Reading in the graph." << endl;
		graph g(vertices, edge_count);
		while (getline(myfile, line)) {
			// Get start, end and weight of an edge.
			int source = atol(line.substr(0, line.find(delimiter)).c_str());
			line = line.substr(line.find(delimiter) + 1, line.size());

			int end = atol(line.substr(0, line.find(delimiter)).c_str());

			int weight = atol(line.substr(line.find(delimiter) + 1, line.size()).c_str());

			g.add(source, end, weight);
		}
		myfile.close();

        // using built-in random generator:
        std::srand ( unsigned ( std::time(0) ) );
        std::random_shuffle ( g.edges.begin(), g.edges.end() );

		int dist[g.total_vertices];
		for (int i=0; i<g.total_vertices;++i) {
			dist[i] = INT_MAX;
		}

		// BELLMAN FORD STARTS HERE.
		dist[source_node] = 0;
        //omp_set_dynamic(0);
        omp_set_num_threads(atoi(argv[2]));
        cout << "Starting bellman ford." << endl;
       	
        double start_time;
        start_time = omp_get_wtime();
        for (int i = 0; i < (g.total_vertices - 1); i++)
		{
            bool any_relaxed = false; 
		    #pragma omp parallel for
			for (int j = 0; j < g.total_edges; j++)
			{
				int u = g.edges[j]->src;
				int v = g.edges[j]->dest;

				int weight = g.edges[j]->weight;

				int prev_dist = dist[v];
                int* mem = &(dist[v]);

                again:
				if (dist[u] != INT_MAX && dist[u] + weight < prev_dist) {
					int new_dist = dist[u] + weight;

//					dist[v] = new_dist;
//					if (0) {
//						goto again;
//					}
                    int old_value =__sync_val_compare_and_swap(mem, prev_dist, new_dist);
                    if (old_value == prev_dist) {
                    	any_relaxed = true;
                    } else {
                    	prev_dist = dist[v];
                    	goto again;
                    }
				}
			}
			// Atleast one node was relaxed in the last loop.
			if (!any_relaxed) {
				break;
			}
		}
		double end_time;
		end_time = omp_get_wtime();
		// BELLMAN FORD ENDS HERE.

		cout << "Time taken for parallel(sec): " << end_time - start_time << endl;

		ofstream myfile;
		myfile.open (filename.c_str());
		for (int i=0; i < g.total_vertices; i++) {
			myfile << i << "," << dist[i] << endl;
		}
		myfile.close();
	}

	return 0;
}


