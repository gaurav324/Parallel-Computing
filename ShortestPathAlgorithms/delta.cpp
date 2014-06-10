/*
 * delta.cpp
 *
 *  Created on: Oct 8, 2013
 *      Author: gnanda
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <set>

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

typedef pair<vector<end*>*, vector<end*>* > heavy_light;

class graph
{
public:
	// Total vertices and edges.
	int total_vertices;
	int total_edges;

	vector< heavy_light* > vertices;

	graph(int vertices, int edge_count) {
		total_vertices = vertices;
		total_edges = edge_count;
	}

	void add(int source, int dest, int weight, int delta) {
		if (source >= vertices.size()) {
			vertices.resize(source + 1);
		}

		if (dest >= vertices.size()) {
			vertices.resize(dest + 1);
		}

		end* e = new end(dest, weight);

		if (vertices[source] == NULL) {
			vertices[source] = new heavy_light(new vector<end*>(), new vector<end*>());
		}

		if (vertices[dest] == NULL) {
			vertices[dest] = new heavy_light(new vector<end*>(), new vector<end*>());
		}

		if (weight > delta) {
			(*(vertices[source])).first->push_back(e);
		} else {
			(*(vertices[source])).second->push_back(e);
		}
	}
};


void relax(vector<set<int>*> *B, vector<int> *dist, int v, int x, int delta) {
	if (x < (*dist)[v]) {

		//cout << 1 << endl;
		int i1 = (*dist)[v] / delta;
		if (i1 <= B->size()) {
			if ((*B)[i1] != NULL) {
				(*B)[i1]->erase(v);
			}
		}

		//cout << 2 << endl;
		int i2 = x / delta;
		//cout << x << "/" << delta << endl;
		if (i2 >= B->size()) {
			//cout << "Resized to: " << i2 + 1 << endl;
			B->resize(i2 + 1, NULL);
		}
		//cout << 3 << endl;
		if ((*B)[i2] == NULL) {
			(*B)[i2] = new set<int>();
		}
		//cout << 4 << endl;
		(*B)[i2]->insert(v);

		//cout << "Relaxing weight of " << v << " to " << x << endl;
		(*dist)[v] = x;
	}
}

int main(int argc, char* argv[]) {

	int delta = atol(argv[2]);
	int s = atol(argv[3]);
	if (delta == 0) {
		cout << "Delta cannot be ZeRo !!" << endl;
		return 0;
	}
	string filename = argv[4];

	cout << "Delta: " << delta << endl;
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

			G.add(source, end, weight, delta);
		}
		myfile.close();

		cout << "Total vertices: " << G.total_vertices << endl;
		cout << "Total edges: " << G.total_edges << endl;

		vector<int> dist;
		for (int i=0; i < G.total_vertices;++i) {
			dist.push_back(INT_MAX);
		}

//		for (int i=0; i < G.total_vertices;++i) {
//			heavy_light* x = G.vertices[i];
//			cout << "Vertex: " << i << endl;
//			for (int j=0; j< x->first->size(); ++j) {
//				cout << (*(x->first))[j]->dest << " " << (*(x->first))[j]->weight << endl;
//			}
//			cout << endl;
//		}

		double start_time;
		start_time = omp_get_wtime();

		cout << "Delta stepping starts here. " << endl;

		vector<set<int>*> B;
		relax(&B, &dist, s, 0, delta);
		int i=0;

//		set<int>* temp = B[0];
//		for (std::set<int>::iterator it=temp->begin(); it!=temp->end(); ++it)
//			cout << ' ' << *it;

//		cout << endl;
		while(i != B.size()) {
			set<int> S;

			//cout << 1 << endl;
			// Do the light ones.
			while(B[i] != NULL && B[i]->size() != 0) {
				vector<pair<int, int> > Req;

				//cout << 2 << endl;
				std::set<int>::iterator it1;
				for(it1 = B[i]->begin(); it1 != B[i]->end(); ++it1) {
					int v = *it1;
					S.insert(v);
                    //cout << "Vertex is: " << v << endl;
					vector<end*>* light = G.vertices[v]->second;
                    //cout << "Light is " << light << endl;
					for(int j=0; j<light->size();j++) {
						int w = (*light)[j]->dest;
						int c = (*light)[j]->weight;

						//cout << 4 << endl;
						Req.push_back(pair<int, int>(w, dist[v] + c));
					}
				}

				//cout << 5 << endl;
				B[i]->clear();

				for (int j=0; j<Req.size(); ++j) {
					pair<int, int> vx = Req[j];
					relax(&B, &dist, vx.first, vx.second, delta);
				}
				//cout << 6 << endl;
			}

			// Do the heavy ones.
			std::set<int>::iterator it1;
			vector<pair<int, int> > Req;
			for(it1 = S.begin(); it1 != S.end(); ++it1) {
				int v = *it1;
				vector<end*>* heavy = G.vertices[v]->first;
				for(int j=0; j < heavy->size();j++) {
					int w = (*heavy)[j]->dest;
					int c = (*heavy)[j]->weight;

					Req.push_back(pair<int, int>(w, dist[v] + c));
				}
			}
			//cout << 7 << endl;

			for (int j=0; j<Req.size(); ++j) {
				pair<int, int> vx = Req[j];
				relax(&B, &dist, vx.first, vx.second, delta);
			}

			//cout << 8 << endl;
			i++;
			//cout << "new I is: " <<  i << endl;
			//cout << "Size of B is: " << B.size() << endl;
		}

		double end_time;
		end_time = omp_get_wtime();

		cout << "Time taken for serial(sec): " << end_time - start_time << endl;

		ofstream myfile;
		myfile.open (filename.c_str());
		for (int i=0; i < G.total_vertices; i++) {
			myfile << i << "," << dist[i] << endl;
		}
		myfile.close();
	}
}
