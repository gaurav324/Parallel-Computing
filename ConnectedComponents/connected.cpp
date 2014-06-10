#include "Galois/Galois.h"
#include "Galois/Graphs/LCGraph.h"
#include "Galois/Statistic.h"
#include "Galois/LargeArray.h"
#include "Galois/Bag.h"
#include "Galois/ParallelSTL/ParallelSTL.h"

#include "Lonestar/BoilerPlate.h"

#include "boost/optional.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <queue>

using namespace std;

static const char* name = "Homework 3: ConnectedComponent ";
static const char* desc = "Computes the Connected Components in the Graph.";
static const char* url = "ConnectedComponents";

enum Algo {
  serial,
  parallel,
};
enum Schedule {
  FIFO,
  ChunkedFIFO,
  dChunkedFIFO,
  OrderedByIntegerMetric,
  LocalQueues,
};

// Parse command line options.
namespace cll = llvm::cl;
static cll::opt<std::string> filename(cll::Positional, cll::desc("<input file>"), cll::Required);
static cll::opt<Algo> algo( cll::desc("Choose an algorithm:"),
                            cll::values(
                                clEnumVal(serial, "Serial"),
                                clEnumVal(parallel, "Parallel"),
                                clEnumValEnd), 
                            cll::init(parallel)
                          );
static cll::opt<Schedule> schedule( cll::desc("Choose a schedule:"),
                                    cll::values(
                                        clEnumVal(FIFO, "FIFO"),
                                        clEnumVal(ChunkedFIFO, "chunked"),
                                        clEnumVal(dChunkedFIFO, "dChunkedFIFO"),
                                        clEnumVal(OrderedByIntegerMetric, "ordered"),
                                        clEnumVal(LocalQueues, "localQ"),
                                        clEnumValEnd), 
                                    cll::init(FIFO)
                                   );

// Create the graph structure.
struct Node {
	int id;
	int compId;
};

typedef Galois::Graph::LC_CSR_Graph<Node, void> Graph;
typedef Graph::GraphNode GNode;

Graph g;

// Magic function used by various Worklists created by Galois system.
struct Process {
    void operator()(GNode n, Galois::UserContext<GNode>& ctx) {
        Node& src_node = g.getData(n, Galois::NONE);
        cout << "Processing: " << src_node.id << endl;
        for(Graph::edge_iterator ii = g.edge_begin(n, Galois::CHECK_CONFLICT), ei = g.edge_end(n, Galois::CHECK_CONFLICT); ii != ei; ++ii) {
            GNode dest = g.getEdgeDst(ii);
            Node& nodeData = g.getData(dest); 

            if (src_node.compId < nodeData.compId) {
                nodeData.compId = src_node.compId;
                ctx.push(dest);
            }
//            // We would try to use compare and swap and then we should be done.
//            int* old_compId_address = &nodeData.compId;
//            while(true) {
//               int prev_value = nodeData.compId;
//               if (src_node.compId < prev_value) {
//                    ctx.push(dest);
//                    if (__sync_bool_compare_and_swap(old_compId_address, prev_value, src_node.compId)) {
//                        break; 
//                    }
//                } else {
//                    break;
//                }
//            }
        }
    }
};

// This is for specifying order in the orderedByIntegerMetric.
struct NodeIndexer
  : std::binary_function<GNode, unsigned int, unsigned int> {
  unsigned int operator() (const GNode& val) const {
    Node& n1_node = g.getData(val);
    return n1_node.compId % 1000;
  }
};

// Main function for executing parallel version of the code.
struct GaloisCCAlgo {
    void operator()() {
        switch(schedule) {
            case FIFO:
                cout << "FIFO worklist enabled." << endl;
                Galois::for_each<GaloisRuntime::WorkList::FIFO<GNode, true> > (g.begin(), g.end(), Process());
                break;
             case ChunkedFIFO:
                cout << "ChunkedFIFO worklist enabled." << endl;
                Galois::for_each<GaloisRuntime::WorkList::ChunkedFIFO<1024> > (g.begin(), g.end(), Process());
                break;
             case dChunkedFIFO:
                cout << "dChunkedFIFO worklist enabled." << endl;
                Galois::for_each<GaloisRuntime::WorkList::dChunkedFIFO<1024> > (g.begin(), g.end(), Process());
                break;
             case OrderedByIntegerMetric:
                cout << "OrderedByIntegerMetric worklist enabled." << endl;
                Galois::for_each<GaloisRuntime::WorkList::OrderedByIntegerMetric<NodeIndexer, GaloisRuntime::WorkList::dChunkedFIFO<1024> > > (g.begin(), g.end(), Process());
                break;
             case LocalQueues:
                cout << "LocalQueues worklist enabled." << endl;
                Galois::for_each<GaloisRuntime::WorkList::LocalQueues<> > (g.begin(), g.end(), Process());
                break;
       }
    }
};

struct SerialFIFO {
    void operator()() {
        queue<GNode> worklist;
        for(auto ii = g.begin(), ie=g.end(); ii != ie; ++ii) {
            worklist.push(*ii);
        }

        while(!worklist.empty()) {
            GNode n = worklist.front();
            worklist.pop();

            auto &src_node = g.getData(n, Galois::NONE);
            for(auto ii = g.edge_begin(n, Galois::NONE), ei = g.edge_end(n, Galois::NONE); ii != ei; ++ii) {
                GNode dest = g.getEdgeDst(ii);
                auto &dest_node = g.getData(dest);
                if (src_node.compId < dest_node.compId) {
                    dest_node.compId = src_node.compId;
                    worklist.push(dest);
                }
            }
        }
    }
};

struct Compare {
	bool operator() (const GNode& n1, const GNode& n2)const {
		Node& n1_node = g.getData(n1);
		Node& n2_node = g.getData(n2);
		return(n1_node.compId > n2_node.compId);
	}
};

struct SerialOrdered {
    void operator()() {
        priority_queue<GNode, vector<GNode>, Compare > worklist;
        for(auto ii = g.begin(), ie=g.end(); ii != ie; ++ii) {
            worklist.push(*ii);
        }

        while(!worklist.empty()) {
            GNode n = worklist.top();
            worklist.pop();

            auto &src_node = g.getData(n, Galois::NONE);
            for(auto ii = g.edge_begin(n, Galois::NONE), ei = g.edge_end(n, Galois::NONE); ii != ei; ++ii) {
                GNode dest = g.getEdgeDst(ii);
                auto &dest_node = g.getData(dest);
                if (src_node.compId < dest_node.compId) {
                    dest_node.compId = src_node.compId;
                    worklist.push(dest);
                }
            }
        }
    }
};

// Switch case to select appropriate version of the serial connected component algorithm.
struct SerialCCAlgo {
    void operator()() {
        switch(schedule) {
            case FIFO:
                cout << "Using Serial-FIFO schedule." << endl;
                SerialFIFO()();
                break;
            case OrderedByIntegerMetric:
                cout << "Using OrderedByComponentId-FIFO schedule. " << endl;
                SerialOrdered()();
                break;
        }
    }
};

int main(int argc, char** argv)
{
    Galois::StatManager statManager;
    LonestarStart(argc, argv, name, desc, url);

    // Initialize timer.
    Galois::StatTimer execTime("ExecutionTime");
    execTime.start();

    // Read the graph into file.
    g.structureFromFile(filename);

    unsigned int id = 0;
    for(Graph::iterator ii = g.begin(), ei = g.end(); ii !=ei; ++ii, ++id) {
        Node& n = g.getData(*ii);
        n.id = id;
        n.compId = id;
    }
    cout << "Total nodes: " << g.size() << endl;

    switch(algo) {
        case serial: 
            SerialCCAlgo()();
            break;
        default: 
            GaloisCCAlgo()(); 
            break;
    }

    for(Graph::iterator ii = g.begin(), ei = g.end(); ii != ei; ++ii) {
        Node& n = g.getData(*ii);
        if (n.compId != 0) {
            //cout << "Non-zero compId for: " << n.id << " and ID is: " << n.compId << endl;
        }
    }

    // Stop timer.
    execTime.stop();

    std::cout << "Time taken to execute: " << execTime.get() << " ms" << endl;
    return 0;
}
