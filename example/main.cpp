#include <iostream>
#include <memgraph/Graph.h>

void help() {
	std::cout << "memgraph file.osm.pbf" << std::endl;
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		help();
		return -1;
	}	
	std::string inFileName( argv[1] );
	memgraph::Graph g( memgraph::Graph::fromPBF(inFileName) );
	
	if (!g.selfCheck()) {
		std::cerr << "Graph is broken" << std::endl;
		return -1;
	}

	std::cout << "Nodes: " << g.nodeCount() << std::endl;
	std::cout << "Edges: " << g.edgeCount() << std::endl;
	return 0;
}
