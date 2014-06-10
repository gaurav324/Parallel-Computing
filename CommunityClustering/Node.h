#include <map>
#include <list>

#include <iostream>

#include "Galois/Galois.h"
#include "Galois/Graph/LCGraph.h"
#include "Galois/Graph/Graph.h"
#include "Galois/Statistic.h"


#include "Galois/Graph/LC_Morph_Graph.h"
#include "Galois/UserContext.h"


#include "Lonestar/BoilerPlate.h"

//#include "Community.h"

////////////////////////////////
using namespace std;

class Cluster_Node
{
// Ignoring scopes for this project.
public:
	Cluster_Node(int _id = 0, int _w_inside_c = 0, int _total_c = 0, int _total_edge_weight = 0) {
		id = _id;
		w_inside_c = _w_inside_c;
		total_c = _total_c;
		total_edge_weight = _total_edge_weight;
	}

    // CLuster_Node id.
    int id;

	// Sum of weights of links inside C.
	int w_inside_c; // Might not be required.

	// Sum of the weights of the links incident to nodes in C.
	int total_c;

	// To form links to the neighboring clusters.
	map<Cluster_Node*, int> cluster_map;

    Galois::Graph::LC_Morph_Graph<Cluster_Node, int>::GraphNode gNode;

	int total_edge_weight;
};

/*******************Graph of Cluster Nodes*************************/
// Graph has Node data, void edge data.
typedef Galois::Graph::LC_Morph_Graph<Cluster_Node, int> Cluster_Graph;

// Opaque pointer to graph node.
typedef Cluster_Graph::GraphNode Cluster_GNode;

class Node
{
// Ignoring scopes for this project.
public:
	Node(int _id = 0, int _loop_weight = 0, int _degree = 0, Cluster_GNode _my_community_GNode = NULL, Cluster_Node* _my_community_node = NULL, int _n_total = 0) {
		id = _id;
		loop_weight = _loop_weight;
		degree = _degree;
		my_community_GNode = _my_community_GNode;
		my_community_node = _my_community_node;
		n_total = _n_total;
	}
    // Node id.
    int id;

    // loop weight.
    int loop_weight;

    // degree of the node, including itself.
    int degree;

    int n_total;

    /***************** community data ..each node is pointing initially to the upper graph*********************/
    //my community
    Cluster_GNode my_community_GNode; // to take locks where necessary
    Cluster_Node* my_community_node;  // to read-only and not to take locks.

	Cluster_GNode leaveExistingCommunity(int node_weight_inside_c, Cluster_Graph& cluster_graph) {

		Cluster_Node& _my_community_node = cluster_graph.getData(my_community_GNode, Galois::MethodFlag::NONE);
		_my_community_node.w_inside_c -= node_weight_inside_c;
		_my_community_node.total_c -= n_total;
		this->my_community_node = NULL;
		//Cluster_GNode temp = my_community_GNode;
		return my_community_GNode;
	}
};

