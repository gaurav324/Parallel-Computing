class Community
{
// Ignoring scopes for this project.
public: 
    // Sum of weights of links inside C.
    int w_inside_c; // Might not be required.

    // Sum of the weights of the links incident to nodes in C.
    int total_c;

    // Constructor.
    Community(int initial_node_loop, int initial_node_degree) {
        // For the stage-1, initial_node_loop would be 1. 
        w_inside_c = initial_node_loop;

        // Initially total_c would just be equal to the
        // degree of the single node.
        total_c = initial_node_degree;
    }
};
