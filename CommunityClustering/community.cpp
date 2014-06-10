#include <iostream>

#include "Galois/Galois.h"
#include "Galois/Graph/LCGraph.h"
#include "Galois/Graph/Graph.h"
#include "Galois/Statistic.h"


#include "Galois/Graph/LC_Morph_Graph.h"
#include "Galois/UserContext.h"


#include "Lonestar/BoilerPlate.h" 

#include "Node.h"

////////////////////
using namespace std;
////////////////////

// Command line options.
namespace cll = llvm::cl;
static cll::opt<std::string> filename(cll::Positional, cll::desc("<input file>"), cll::Required);

// Constants for the Galois program.
static const char* name = "Machine learing - Community Clustering. ";
static const char* desc = "Finds Community clusters.";
static const char* url = "Community Clusters.";

/********************Graph of Nodes ********************************/
// Graph has Node data, void edge data.
typedef Galois::Graph::LC_Morph_Graph<Node, int> Graph;

// Opaque pointer to graph node.
typedef Graph::GraphNode GNode;


// Sum total of all the nodes in the graph.
int total_weight = 0;

/*******************************************************************/
// This would try to find the local maxima for the given set of nodes.
struct levelOneProcess {
	Graph& g;
	Cluster_Graph& cluster_g;
	levelOneProcess(Graph& _g, Cluster_Graph& _cluster_g):g(_g),cluster_g(_cluster_g) {}
    void operator()(GNode n, Galois::UserContext<GNode>& ctx) {
    	Node& src_node = g.getData(n, Galois::MethodFlag::NONE);
        Cluster_GNode src_Gcomm = src_node.my_community_GNode;
		Cluster_Node* src_comm = &cluster_g.getData(src_Gcomm, Galois::MethodFlag::CHECK_CONFLICT);

		// A map to hold a sum of weight of links from i inside c.
        std::map<Cluster_Node*, int, less<Cluster_Node*>, Galois::PerIterAllocTy::rebind<pair<const Cluster_Node*, int> >::other> node_weight_inside_c(less<Cluster_Node*>(), ctx.getPerIterAlloc());

		//map<Cluster_Node*, int> node_weight_inside_c; //per iteration allocation gmetis..
		//map<Cluster_Node*, Cluster_GNode> mp;
        //std::map<Cluster_Node*, Cluster_GNode, less<Cluster_Node*>, Galois::PerIterAllocTy::rebind<std::pair<Cluster_Node*,Cluster_GNode> >::other> mp(less<Cluster_Node*>(), ctx.getPerIterAlloc());

        //mp[src_comm] = src_Gcomm;
        node_weight_inside_c[src_comm] = 0;
		for (auto ii = g.edge_begin(n, Galois::MethodFlag::NONE),ei = g.edge_end(n, Galois::MethodFlag::NONE); ii != ei; ++ii) {
			GNode dest = g.getEdgeDst(ii);
			Node& temp = g.getData(dest, Galois::MethodFlag::NONE);

			Cluster_Node* neigh_community = &cluster_g.getData(temp.my_community_GNode, Galois::MethodFlag::CHECK_CONFLICT);
			node_weight_inside_c[neigh_community] +=  g.getEdgeData(ii); // 2 links added a---b... one for a and one for b
		}

		// Remove the node in context from the community to which it belongs and see if we get gain.
		Cluster_GNode leftComm = src_node.leaveExistingCommunity(node_weight_inside_c[src_comm] * 2 + src_node.loop_weight, cluster_g);

		// For all the adjacent nodes of the n, find which one has the profit of modularity.
		double maxProfit = 0;
		Cluster_GNode maxProfitGComm = leftComm;

		for(auto iterator = node_weight_inside_c.begin(); iterator != node_weight_inside_c.end(); ++iterator) {
			Cluster_Node* neigh = iterator->first;
			double new_profit = iterator->second  - (neigh->total_c * src_node.n_total * 1.0 / total_weight);
			if (new_profit > maxProfit) {
				maxProfit = new_profit;
				//maxProfitGComm = mp[neigh]; // start pointing to the neighbor's community.
				maxProfitGComm = neigh->gNode; // start pointing to the neighbor's community.
				if (neigh->gNode == NULL) {
					cout << "Assigned " << neigh->gNode << " NULL value." << endl;
				}
			}
		}

		if (maxProfitGComm != NULL) {
            //cout << maxProfitGComm	 << endl;
		} else {
			cout << "NULL maxProfit" << endl;
		}
		Cluster_Node* maxProfitComm = &cluster_g.getData(maxProfitGComm, Galois::MethodFlag::NONE);
		maxProfitComm->w_inside_c += node_weight_inside_c[maxProfitComm] * 2 + src_node.loop_weight;
		maxProfitComm->total_c += src_node.n_total;
        src_node.my_community_node = maxProfitComm;
		src_node.my_community_GNode = maxProfitGComm;
	}
};

/****** To calculate loop_weight and edge weights to other communities ******/
struct PhaseTwo_pre {
	Graph& graph;
	Cluster_Graph& cluster_g;
	PhaseTwo_pre(Graph& _g, Cluster_Graph& _cluster_g):graph(_g),cluster_g(_cluster_g){}

	void operator()() {
	for(auto ii = graph.begin(); ii != graph.end(); ii++) {
		Node& src_node = graph.getData(*ii);

		for(auto jj = graph.edge_begin(*ii); jj != graph.edge_end(*ii); ++jj) {
			GNode node = graph.getEdgeDst(jj);
			Node& dest_node = graph.getData(node);

			if(src_node.my_community_node != dest_node.my_community_node) {
				src_node.my_community_node->cluster_map[dest_node.my_community_node] +=  graph.getEdgeData(jj);
				src_node.my_community_node->total_edge_weight +=  graph.getEdgeData(jj);
			}
		}
	  }
	}
};

double modularity(Graph& g, Cluster_Graph& cluster_g) {
    double mod = 0;

    map<Cluster_Node*, int> comm;
    for(Graph::iterator ii = g.begin(), ei = g.end(); ii != ei; ++ii) { 
        Node& n = g.getData(*ii);
        comm[n.my_community_node] = 1;
    }

    // cout << "**********************Size:" << comm.size() << endl;
    for(auto iterator = comm.begin(); iterator != comm.end(); ++iterator) {
    	Cluster_Node* neigh = iterator->first;
    	//Cluster_Node& neigh = cluster_g.getData(neigh);

    	if (neigh->total_c > 0) {
    		mod += (neigh->w_inside_c * (double)1.0/total_weight) - ((neigh->total_c * (double)1.0/total_weight) * (neigh->total_c * (double)1.0/total_weight));
    	}
    }   
    //cout << "Modularity is: " << mod << endl;
    return mod;
}

double performLevelOne(Graph& g, Cluster_Graph& cluster_g){
        double old_mod = modularity(g, cluster_g);
        double new_mod = old_mod;

        //cout << "Inside performLevelOne \n";
        for (auto ii = g.begin(); ii != g.end(); ++ii) {
           Node& no = g.getData(*ii); 
        }
        do {
            old_mod = new_mod;
            typedef Galois::WorkList::ChunkedFIFO<1024> WL; //worklist to be used.
            Galois::for_each(g.begin(), g.end(), levelOneProcess(g, cluster_g), Galois::wl<WL>());
            new_mod = modularity(g, cluster_g);
            //cout << "Internal to level one. New mod: " << new_mod << endl;
        }
        while((new_mod - old_mod) > 0.00001);

        // Preprocessing for new graph construction.
        PhaseTwo_pre(g, cluster_g)();
        return new_mod;
 }

void makeNewGraphs(Cluster_Graph& cluster_graph_new, Graph& graph_new, Cluster_Graph& cluster_g) {
		//cout << "Building new graphs." << endl;
	    int id = 0;
	    map<Cluster_Node*, Cluster_GNode> map_for_edges_cluster;
	    map<Cluster_Node*, GNode> map_for_edges_graph;

        //int weight = 0;
        //int supposed_edges = 0;
	    for(auto ii = cluster_g.begin(); ii != cluster_g.end(); ii++) {
	    	Cluster_Node& c_node = cluster_g.getData(*ii);
	    	if (c_node.total_c > 0) {
	    		int numEdges = c_node.cluster_map.size();
                //supposed_edges += numEdges;
	    		Cluster_GNode cluster_GNode = cluster_graph_new.createNode(numEdges, id, c_node.w_inside_c, c_node.total_edge_weight + c_node.w_inside_c, 0);
                map_for_edges_cluster[&c_node] = cluster_GNode;
	    		Cluster_Node& my_community = cluster_graph_new.getData(cluster_GNode);
                my_community.gNode = cluster_GNode;
	    		GNode Graph_N = graph_new.createNode(numEdges, id, c_node.w_inside_c, numEdges + c_node.w_inside_c, 
                                                     cluster_GNode, &my_community, c_node.total_edge_weight + c_node.w_inside_c);
                //weight += c_node.w_inside_c;
	    		map_for_edges_graph[&c_node] = Graph_N ;
	    		id++;
            }
	    }
	    //cout << "map_for_edges_cluster = " << map_for_edges_cluster.size() << " map_for_edges_graph = " << map_for_edges_graph.size() << "\n";

       //int edges = 0;
       //int edge_weights = 0;
	   for(auto ii = cluster_g.begin(); ii != cluster_g.end(); ii++) {
	    	Cluster_Node& c_node = cluster_g.getData(*ii);
	    	if(c_node.total_c > 0) {
	    		//int c_node_id = c_node.id;

	    		Cluster_GNode& src_cluster = map_for_edges_cluster[&c_node];
	    		GNode& src_graph = map_for_edges_graph[&c_node];

	    		auto make_edge_ii = c_node.cluster_map.begin();
	    		for(; make_edge_ii != c_node.cluster_map.end(); make_edge_ii++) {
	    			Cluster_Node* dest_node = (*make_edge_ii).first;
	    			//int dest_id = (*dest_node).id;
	    			int edge_weight = (*make_edge_ii).second;

	    			Cluster_GNode& dest_cluster = map_for_edges_cluster[dest_node];
                    cluster_graph_new.addEdge(src_cluster, dest_cluster, Galois::MethodFlag::NONE, edge_weight);

	    			GNode& dest_graph = map_for_edges_graph[dest_node];
	    			graph_new.addEdge(src_graph, dest_graph, Galois::MethodFlag::NONE, edge_weight);
                    //weight += edge_weight;
                    //edges += 1;
                    //edge_weights += edge_weight;
	    		}
	    	}
	    }
        //cout << "Supposed Edges: " << supposed_edges << endl;
        //cout << "Total edges: " << edges << endl;
        //cout << "Total weight: " << weight << endl;
        //cout << "Total edge weights: " << edge_weights << endl;
}

//// START OF THE MAIN FUNCTION ////
int main(int argc, char** argv) {
    Galois::StatManager statManager;
    LonestarStart(argc, argv, name, desc, url);

    cout << "******starting*********" << endl;

    // Instantiate the graph.
    Graph* Node_graph = new Graph();
    Cluster_Graph* Cluster_graph = new Cluster_Graph(); // another graph for clusters. so that we can delete.

    // Read graphs from the file.
    Galois::Graph::readGraph(*Node_graph, filename);
    Galois::Graph::readGraph(*Cluster_graph, filename);

    //Making Cluster Graph
    int id = 0;
    for(auto ii = (*Cluster_graph).begin(); ii != (*Cluster_graph).end(); ++ii, ++id) {
    	Cluster_Node& n  = (*Cluster_graph).getData(*ii);

    	n.id = id;
    	n.w_inside_c = 0;
    	n.total_edge_weight = 0;
        n.gNode = *ii;

    	// Intializing edge weights.
    	for(auto jj = (*Cluster_graph).edge_begin(*ii) , kk = (*Cluster_graph).edge_end(*ii); jj != kk; jj++) {
    	    (*Cluster_graph).getEdgeData(jj) = 1;
    		n.total_c += (*Cluster_graph).getEdgeData(jj); // It is the total weight incident on a node in a cluster.
    	}
    }

    cout << "Nodes: = " << id << endl;

    /******************************* Making Node graph ******************************/
    id = 0;
    auto c_ii = (*Cluster_graph).begin(); // to map lower graph to the upper cluster graph

    for(auto ii = (*Node_graph).begin(), ei = (*Node_graph).end(); ii != ei; ++ii, ++id, ++c_ii) {
        Node& n = (*Node_graph).getData(*ii);
        Cluster_Node& c_n = (*Cluster_graph).getData(*c_ii);

        n.id = id;
        n.loop_weight = 0;
        
        // This would include the self loop.
        n.degree = (int)std::distance((*Node_graph).edge_begin(*ii, Galois::NONE), (*Node_graph).edge_end(*ii, Galois::NONE));

        const Cluster_GNode Gnode_ref = *c_ii;
        n.my_community_GNode = Gnode_ref;
        n.my_community_node = &c_n;

        //Intializing edge weights.
		 for(auto jj = (*Node_graph).edge_begin(*ii),kk=(*Node_graph).edge_end(*ii);jj!=kk;jj++) {
			 (*Node_graph).getEdgeData(jj)=1;
			 total_weight += (*Node_graph).getEdgeData(jj);
			 n.n_total += (*Node_graph).getEdgeData(jj);
			  // weight+=1;
			}
    }
	cout << "Total weight = " << total_weight << "\n";

    /****************************************************************/
    // Perform the first level on the graph, which would 
    // produce clusters.

    Cluster_Graph* Current_Cluster_graph = Cluster_graph;
    Graph* Current_Node_graph = Node_graph;
    Graph* prev_Node_graph = Node_graph;
    Cluster_Graph* prev_cluster_graph = Cluster_graph;
    double prev_mod = modularity(*Current_Node_graph,*Current_Cluster_graph);
    double new_mod = prev_mod;

    // Initialize timer.
    Galois::StatTimer execTime("ExecutionTime");
    execTime.start();

    cout << "Initially modularity is: " << new_mod << endl;
    do {
    		prev_mod = new_mod;
    	    new_mod = performLevelOne(*Current_Node_graph, *Current_Cluster_graph);
    	    cout << "Modularity after phase one calculation = " << new_mod << endl;
    	    if (new_mod <= prev_mod) {
    	    	cout << endl << "End of the iterations." << endl;
    	    	break;
    	    }
    	    Cluster_Graph* Cluster_graph_new = new Cluster_Graph();
			Graph* graph_new = new Graph();

			prev_Node_graph = Current_Node_graph;
			prev_cluster_graph = Current_Cluster_graph;
            for (auto ii = (*prev_Node_graph).begin(); ii != (*prev_Node_graph).end(); ++ii) {
                Node& node = (*prev_Node_graph).getData(*ii);
                if (node.my_community_node == NULL) {
                    cout << "Node has no community." << endl; 
                }
            }
            int total_weight = 0;
            int total_w_inside = 0;
            /*
             * *for (auto ii = (*prev_cluster_graph).begin(); ii != (*prev_cluster_graph).end(); ++ii) {
                Cluster_Node& node = (*prev_cluster_graph).getData(*ii);
                //if (node.w_inside_c < 0) {
                //    cout << "Weight is ." << node.w_inside_c << endl; 
                //}
                total_weight +=  node.total_c;
                total_w_inside += node.w_inside_c;
            }
            cout << "Sum of weights inside is: " << total_w_inside << endl;
            cout << "Total weight of cluster finally is: " << total_weight << endl;
			*/
            delete prev_Node_graph;
			makeNewGraphs(*Cluster_graph_new, *graph_new, *Current_Cluster_graph);
			delete prev_cluster_graph;
			int nodes = (int)distance((*Cluster_graph_new).begin(), (*Cluster_graph_new).end());
			cout << endl << "New number of nodes = " << nodes << endl;
			Current_Node_graph = graph_new;
			Current_Cluster_graph = Cluster_graph_new;

			//cout << "Calling performLevelOne again." << endl;

			//performLevelOne(graph_new);
    //} while(0);
    } while((new_mod - prev_mod) > 0.00001);

    // Stop timer.
    execTime.stop();

    std::cout << "Time taken to execute: " << execTime.get() << " ms" << endl;

    /*for(auto ii = Cluster_g.begin();ii != Cluster_g.end(); ++ii) {
    	Cluster_Node c_node = Cluster_g.getData(*ii);
    	if(c_node.w_inside_c!=0){
    		cout<<"ID = "<<c_node.id<<"\t loop_weight = "<<c_node.w_inside_c<<endl;
    	}
    }*/

   /******************checking new graphs***************************/
/*   for(auto ii = (*graph_new).begin(); ii!=(*graph_new).end(); ++ii) {
	   Node& n  = (*graph_new).getData(*ii);
	   cout<<"new ID = "<<n.id<<"\n\t";
	   for(auto jj = (*graph_new).edge_begin(*ii); jj!= (*graph_new).edge_end(*ii); jj++) {
		   GNode neigh = (*graph_new).getEdgeDst(jj);
		   Node& neigh_data = (*graph_new).getData(neigh);
		   int egdeW = (*graph_new).getEdgeData(jj);
		   cout<<neigh_data.id<<" edgeW = "<<egdeW<<"\t";
	   }
	   cout<<"\n";

   }
*/

/*
   for(auto ii = Cluster_graph_new.begin(); ii!=Cluster_graph_new.end(); ++ii) {
   	   Cluster_Node& n  = Cluster_graph_new.getData(*ii);
   	   cout<<"new ID = "<<n.id<<"\n\t";
   	   for(auto jj = Cluster_graph_new.edge_begin(*ii); jj!= Cluster_graph_new.edge_end(*ii); jj++) {
   		   Cluster_GNode neigh = Cluster_graph_new.getEdgeDst(jj);
   		   Cluster_Node& neigh_data = Cluster_graph_new.getData(neigh);
   		   int egdeW = Cluster_graph_new.getEdgeData(jj);
   		   cout<<neigh_data.id<<" edgeW = "<<egdeW<<"\t";
   	   }
   	   cout<<"\n";

      }

*/
    //printCommunities();
    
    
  //  cout << "***************" << endl;
    return(0);
}
