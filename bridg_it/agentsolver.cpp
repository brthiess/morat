
#include <cmath>
#include <string>
#include <vector>
#include <stack>  

#include "../lib/alarm.h"
#include "../lib/fileio.h"
#include "../lib/string.h"
#include "../lib/time.h"

#include "agentsolver.h"
#include "board.h"


/***************************
 * 
 * This is how the Bridg-It board is translated into a graph
 * Note:  First and last rows are condensed to single point 
 * 
 * 
 * 	    1		First Row
 * 	  / | \
 * 	 2--3--4	Second Row
 * 	 |  |  |
 * 	 5--6--7	Third Row
 * 	  \ | /
 * 	    8		Last Row 
 * 
 ****************************/

namespace Morat {
namespace BridgIt {



/**
 * This function is called when genmove is called
 */
void AgentSolver::search(double time, uint64_t max_runs, int verbose){
	//Get the adjacency list representation of the board
	AgentSolver::Adjacency_List board_matrix = AgentSolver::getAdjacencyList();
	printf("\nBoard Matrix");
	board_matrix.graph_to_s();
	
	//Get rid of problem vertices
	//board_matrix = remove_problem_vertices(board_matrix);
	
	//Find 2 Edge Disjoint Spanning Trees
	std::vector<AgentSolver::Adjacency_List> trees = find_edge_disjoint_trees(board_matrix);
	
	AgentSolver::get_best_move(trees);
}

/**
 * Is given a vector of edge disjoint trees
 * And determines the best move to play
 */
void AgentSolver::get_best_move(std::vector<AgentSolver::Adjacency_List> trees) {
	AgentSolver::Adjacency_List tree1;
	AgentSolver::Adjacency_List tree2;
	
	
	int best_move = -1;
	
	pool.pause();
	
	//Found Win
	if (trees.size() >= 2){
		logerr("Solved as a win\n");

		//Here for testing purposes
		logerr("PV:\n");
		//Grab the trees
		tree1 = trees.at(0);
		tree2 = trees.at(1);		
	}
	else {
		logerr("Solved as a loss\n");
		logerr("PV:\n");
		root.outcome == Outcome::UNKNOWN;
		return;
	}
	
	pool.reset();
	
	//Delete the edge that connects the starting vertice and the finishing vertice
	tree1.delete_edge(0, tree1.get_number_of_vertices() - 1);
	tree2.delete_edge(0, tree2.get_number_of_vertices() - 1);
	

	//Swap trees so that tree2 is not connected and tree1 is
	if (!is_connected(tree1) && is_connected(tree2)){
		printf("\nTree 1 is not connected and tree 2 is");
		Adjacency_List temp_tree;
		temp_tree = tree1;
		tree1 = tree2;
		tree2 = temp_tree;	
	}
	
	//Grab all edges from tree 1
	std::vector<AgentSolver::Edge> edges = get_all_edges(tree1);
	for (int e = 0; (unsigned) e < edges.size(); e++) {
		int v1 = edges.at(e).getV1(); 
		int v2 = edges.at(e).getV2();
		
		if (!vertices_are_in_the_same_set(tree2, v1, v2) && tree1.get_number_of_duplicate_edges(v1, v2) == 1) {
			std::cout << "\n V1 : " << v1 << " and V2: " << v2 << "are not in the same set";
			std::cout << "\nPlay Move: "<< v1 << ", " << v2;
			v1 = tree1.get_original_id(v1);
			v2 = tree1.get_original_id(v2);
			best_move = AgentSolver::edge_to_xy(v1, v2);
			break;
		}
	}
	Move m;
	
	if (best_move != -1) {
		m = rootboard.yx(best_move);
	}
	//Board position is so good that we can win regardless if we play first or second
	else {
		std::cout <<"\n\nPlaying Random Move!";
		m = get_random_move();
		best_move = rootboard.xy(m);
	}
	
	root.outcome = rootboard.toplay();
	root.bestmove = m;	
	
	std::cout << "\n\nBest Move: " << best_move;
}

/**
 * Returns a random move
 */
 
Move AgentSolver::get_random_move() {
	for (int x = 1; x < rootboard.vecsize(); x++) {
		for (int y = 1; y < rootboard.vecsize(); y++) {
			if (rootboard.valid_move(x,y)) {
				return Move(x,y);
			}
		}
	}
	return Move(1,1);
}

/**
 * Will return two edge disjoint trees, unless none are found in which case, 
 * it returns none
 */
std::vector<AgentSolver::Adjacency_List> AgentSolver::find_edge_disjoint_trees(Adjacency_List board_matrix) {
	
	//Create a vector to contain the edge disjoint trees
	std::vector<AgentSolver::Adjacency_List> edge_disjoint_trees;
	
		//1. Create two identical spanning trees (tree1 and tree2)
	//Find a spanning tree from the main graph
	Adjacency_List baseTree = get_spanning_tree(board_matrix);
	
	//Create two identical trees from the base tree
	Adjacency_List tree1 = baseTree;
	Adjacency_List tree2 = baseTree;
	
		//2. Create a graph from remaining edges (called common_chords ?)
	//Create a second tree by removing the edges of baseTree from the main graph
	Adjacency_List leftovers = subtract_trees(board_matrix, baseTree);
	
	//3. For each edge c in common_chords
			//a. add c to tree1
			//b. Find all the edges in the cycle it created
			//c. For each edge r in the cycle
				//i. If r is in both tree1 and tree2
					//1. Delete r from tree1
					//2. Break from loop 
				//ii. If no r's are in both trees
					//1. Repeat Step 3 but add c to t2 instead and delete c from tree1
					
			//d. Delete c from common_chords
	//4. Remove edges that are in both tree1 and tree2 from tree1 and tree2
	
	//5. Remove vertices in tree1 and tree2 that are not connected to their start and end vertices
	
		//6. Add tree1 and tree2 to a vector and return it
	//add edge disjoint trees to edge_disjoint_trees
	edge_disjoint_trees.push_back(tree1);
	edge_disjoint_trees.push_back(tree2);
	
	//return the edge disjoint trees
	return edge_disjoint_trees;
	
}


/**
 * Will return two edge disjoint trees, unless none are found in which case, 
 * it returns only one.
 */
std::vector<AgentSolver::Adjacency_List> AgentSolver::find_edge_disjoint_trees_old(Adjacency_List board_matrix) {	
	//Vector containing all edge disjoint trees found
	std::vector<AgentSolver::Adjacency_List> edge_disjoint_trees;
	
	//Find a spanning tree from the main graph
	Adjacency_List tree1 = get_spanning_tree(board_matrix);
	
	//Subtract the edges the edges from the main graph by tree1 and get a second tree
	Adjacency_List tree2 = subtract_trees(board_matrix, tree1);
	
	//To keep track of the edges used for the swapping
	std::vector<AgentSolver::Edge> edges_used;
	std::vector<AgentSolver::Edge> tree1_edges = AgentSolver::get_all_edges(tree1);
	std::vector<AgentSolver::Edge> tree2_edges = AgentSolver::get_all_edges(tree2);
	
	printf("\n\nOriginal Spanning Tree");
	tree1.graph_to_s();
	
	printf("\n\nSubtracted Tree");
	tree2.graph_to_s();

	while(not_all_edges_used(edges_used, board_matrix)) {
		
		printf("\n\nTree 1");
		tree1.graph_to_s();
	
		printf("\n\nTree 2");
		tree2.graph_to_s();
		//Check if both trees are connected
		//If so then we have found 2 edge disjoint trees
		if (is_connected(tree1) && is_connected(tree2)) {
			edge_disjoint_trees.push_back(tree1);
			edge_disjoint_trees.push_back(tree2);
			printf("\nBoth Are connected!");
			printf("\nConnected Tree 1");			
			tree1.graph_to_s();
			printf("\nConnected Tree 2");
			tree2.graph_to_s();
			break;
		}
		//Else find a connecting edge for tree2 (so that it is fully connected)
		//and swap it with tree1
		else {
			printf("\n\nAt least one tree not connected");
			swap_edges(&tree1, &tree2, &edges_used, &tree1_edges, &tree2_edges);			
		}
	}
	
	printf("\n\nTree 1");
	tree1.graph_to_s();
	printf("\n\nTree 2");
	tree2.graph_to_s();
	return edge_disjoint_trees;	
}

/**
 * Is given a set of edges and vertices (represented by an adjacency list)
 * and returns an adjacency list of a spanning tree
 * Returns null (?) if none exist
 */
AgentSolver::Adjacency_List AgentSolver::get_spanning_tree(Adjacency_List board_matrix) {
	
	//List of all visited vertices
	std::vector<int> visited_vertices;
	
	int number_of_vertices = board_matrix.get_number_of_vertices();
	AgentSolver::Adjacency_List spanning_tree (number_of_vertices);
	
	//Used to keep track of where we are in graph
	std::stack<int> s;
	
	//Push first vertice onto stack
	s.push(0);
	visited_vertices.push_back(0);
	while(!all_vertices_visited(visited_vertices, board_matrix)) {
		while(!s.empty()) {
			int v1 = s.top();
			//Go through all vertices, 
			//check if it is connected to the current vertice 
			//and if it has not been visited before
			for(int v2 = 0; v2 < number_of_vertices; v2++) {
				spanning_tree.set_original_id(v2, board_matrix.get_original_id(v2));
				if (board_matrix.is_connected(v1, v2) && std::find(visited_vertices.begin(), visited_vertices.end(), v2) == visited_vertices.end()) {
					spanning_tree.addEdge(v1,v2);
					spanning_tree.addEdge(v2,v1);
					s.push(v2);
					visited_vertices.push_back(v2);
					v1 = s.top();
					v2 = 0;
				}
			}
			s.pop();
		}
		//Clear the stack
		clear(s);
		

		//Pop the unvisited vertice onto the stack
		for(int v = 0; v < number_of_vertices; v++) {
			//If v is not in visited_vertices array, then it has not been visited
			//Push it onto stack and start depth first search again
			//If we find a vertice here, then we know the graph is not connected
			if (std::find(visited_vertices.begin(), visited_vertices.end(), v) == visited_vertices.end()) {
				s.push(v);
				visited_vertices.push_back(v);
				break;
			}
		}
	}
	return spanning_tree;
}

void AgentSolver::clear( std::stack<int> &s ) {
   std::stack<int> empty;
   std::swap( s, empty );
}


bool AgentSolver::all_vertices_visited(std::vector<int> visited_vertices, Adjacency_List tree) {
	int number_of_vertices = tree.get_number_of_vertices();
	for (int v = 0; v < number_of_vertices; v++) {
		//If a vertice is not found in the visited_vertices array, return false
		if ( std::find(visited_vertices.begin(), visited_vertices.end(), v) == visited_vertices.end()) {
			return false;
		}
	}
	return true;
}

/**
 * Is given 2 trees (represented as adjacency lists).  
 * The second tree's edges are subtracted from the first
 * It returns the result
 */
AgentSolver::Adjacency_List AgentSolver::subtract_trees(Adjacency_List main_tree, Adjacency_List subtracting_tree) {
	
	//Make a copy
	int number_of_vertices = main_tree.get_number_of_vertices();
	AgentSolver::Adjacency_List copy_tree = copyTree(main_tree);

	
	printf("\n***Copy Tree***");
	copy_tree.graph_to_s();
	
	//Go through all possible edge combinations and check to see if both graphs have them
	for (int v1 = 0; v1 < number_of_vertices; v1++) {
		for (int v2 = 0; v2 < number_of_vertices; v2++) {
			//If the main graph is connected and the second graph is not, then add it to the subtracted tree
			if (subtracting_tree.is_connected(v1,v2)) {
				for(int i = 0; i < subtracting_tree.get_number_of_duplicate_edges(v1,v2); i++) {
					//printf("\n\n***Copy Tree Before Deleting: V1: %d V2: %d", v1, v2);
					//copy_tree.graph_to_s();
					copy_tree.delete_edge(v1,v2);
					//printf("\n\n***Copy Tree After Deleting: V1: %d V2: %d", v1, v2);
					//copy_tree.graph_to_s();
				}
			}
		}	
	}
	return copy_tree;
}

/**
 * Checks if the given tree connects every vertice
 * Returns true if it is connected
 */
bool AgentSolver::is_connected(Adjacency_List tree) {
	//List of all visited vertices
	std::vector<int> visited_vertices;
	
	int number_of_vertices = tree.get_number_of_vertices();
	AgentSolver::Adjacency_List spanning_tree (number_of_vertices);
	
	//Used to keep track of where we are in graph
	std::stack<int> s;
	
	//Push first vertice onto stack
	s.push(0);
	visited_vertices.push_back(0);
	
		while(!s.empty()) {
			int v1 = s.top();
			//Go through all vertices, 
			//check if it is connected to the current vertice 
			//and if it has not been visited before
			for(int v2 = 0; v2 < number_of_vertices; v2++) {
				if (tree.is_connected(v1, v2) && std::find(visited_vertices.begin(), visited_vertices.end(), v2) == visited_vertices.end()) {
					spanning_tree.addEdge(v1,v2);
					spanning_tree.addEdge(v2,v1);
					s.push(v2);
					visited_vertices.push_back(v2);
					v1 = s.top();
					v2 = 0;
				}
			}
			s.pop();
		}
		
	//Check to see if we have visited all the vertices through depth first search
	//If we have not, then the graph is not connected
	if (visited_vertices.size() < (unsigned)number_of_vertices) {
		printf("\nTree is not connected");
		return false;
	}
	else {
		printf("\nTree is connected");
		return true;
	}
}

/**
 * Is given two trees (all in Adjacency_List pointer form) 
 * and a vector of already used edges
 * It finds an edge that will connect the unconnected tree
 * and that has not been swapped yet 
 * and swaps the two
 */
void AgentSolver::swap_edges(Adjacency_List *tree1, Adjacency_List *tree2, std::vector<AgentSolver::Edge> *used_edges,
								std::vector<Edge> *tree1_edges, std::vector<Edge> *tree2_edges) {
	Adjacency_List * temp_tree;
	std::vector<Edge> * temp_edges;
	//int number_of_vertices = AgentSolver::get_number_of_vertices();	
	
	//Swap trees so that tree2 is not connected and tree1 is
	if (!is_connected(*tree1) && is_connected(*tree2)){
		printf("\nTree 1 is not connected and tree 2 is");
		temp_tree = tree1;
		tree1 = tree2;
		tree2 = temp_tree;
		temp_edges = tree1_edges;
		tree1_edges = tree2_edges;
		tree2_edges = temp_edges;		
	}
		
	printf("\n\n***Trees Before Swap***");
	printf("\n\nTree 1");
	tree1->graph_to_s();	
	printf("\n\nTree 2");
	tree2->graph_to_s();
		
		//Iterate through all edges in tree 1.
		int e = 0;
		for (; (unsigned) e < tree1_edges->size(); e++) {
			int v1 = tree1_edges->at(e).getV1(); 
			int v2 = tree1_edges->at(e).getV2();
			long id = tree1_edges->at(e).getID();
			printf("\n*****V1: %d V2: %d  ID: %ld", v1, v2, id);
			//Check to make sure edge has not been swapped already
			if (!edge_in(*used_edges, id)) {
				printf("\nEdge is not in used Edges");
				//Check to see if this edge would connect the other tree (tree 2)
				if (!vertices_are_in_the_same_set(*tree2, v1, v2)) {
					printf("Vertices are not in the same set in tree 2");
					//Swap edges
					tree1->delete_edge(v1, v2);
					tree1->delete_edge(v2, v1);
					tree2->addEdge(v1, v2);
					tree2->addEdge(v2, v1);
					used_edges->push_back(tree1_edges->at(e));
					tree2_edges->push_back(tree1_edges->at(e));
					tree1_edges->erase(tree1_edges->begin() + e);
					break;
				}
			}
		}
		//If we iterated over all edges and could not find a suitable edge
		//Pick a random edge instead
		if ( (unsigned) e == tree1_edges->size()) {
			printf("\nNo Edge Found.  Picking Random edge...");

			
			for (int i = 0; (unsigned) i < tree1_edges->size(); i++) {
				int v1 = tree1_edges->at(i).getV1(); 
				int v2 = tree1_edges->at(i).getV2();
				long id = tree1_edges->at(i).getID();
				printf("\n*****Random Vertices: V1: %d V2: %d  ID: %ld", v1, v2, id);
				//Check to see if this edge would connect the other tree (tree 2)
				if (!vertices_are_in_the_same_set(*tree2, v1, v2)) {
					printf("Vertices are not in the same set in tree 2");
					//Swap edges
					tree1->delete_edge(v1, v2);
					tree1->delete_edge(v2, v1);
					tree2->addEdge(v1, v2);
					tree2->addEdge(v2, v1);
					used_edges->push_back(tree1_edges->at(i));
					tree2_edges->push_back(tree1_edges->at(i));
					tree1_edges->erase(tree1_edges->begin() + i);
					break;
				}
			}			
		}
			
	printf("\n\n***Trees After Swap***");
	printf("\n\nTree 1");
	tree1->graph_to_s();
	
	printf("\n\nTree 2");
	tree2->graph_to_s();
			//	1. Check to see if (v1,v2) is an edge in tree1 and it has not already been swapped previously
			//		a. If so,Check to see if (v1,v2) are in the same set for tree2
			//			a. If not, then swap the edges between the trees
			//			b. If they are, then go back to 1. and repeat.
			//		b. If not, then repeat 1.
			
	printf("\n\n\n\n\n\n\n");

}


/**
 * Is given an adjacency matrix and removes the vertices with 
 * 0 or 1 connections to other vertices
 */
 AgentSolver::Adjacency_List AgentSolver::remove_problem_vertices(Adjacency_List tree) {
	 int number_of_vertices = tree.get_number_of_vertices();
	 std::vector<Edge> edges;
	 for (int v = number_of_vertices - 1 ; v >= 0; v--) {
		 if (tree.get_number_of_edges_attached_to_vertice(v) <= 1) {
			 std::cout<<"\nFound a problem vertice at vertice: " << v;
			 tree.delete_vertex(v);
			 v = tree.get_number_of_vertices();
		 }
	 }
	 return tree;
 }

/**
 * Grabs all the edges from an adjacency list
 */ 
std::vector<AgentSolver::Edge> AgentSolver::get_all_edges(AgentSolver::Adjacency_List tree) {
	
	//Grab the number of vertices
	int number_of_vertices = tree.get_number_of_vertices();
	
	//Make a copy
	Adjacency_List copy_tree = copyTree(tree);	

	std::vector<AgentSolver::Edge> edges;
	
	for (int u = 0; u < number_of_vertices; u++) {
		for( int v = u; v < number_of_vertices; v++) {
			if (copy_tree.delete_edge(u,v)) {
				edges.push_back(AgentSolver::Edge(u,v, get_id()));
				v -= 1;
			}
		}
	}

	return edges;
}

/**
 * Returns a copy of the given tree
 */
 
 AgentSolver::Adjacency_List AgentSolver::copyTree(Adjacency_List tree) {
	//Get the number of vertices of the graph
	int number_of_vertices = tree.get_number_of_vertices();
	Adjacency_List copy_tree  (number_of_vertices);
	for (int v1 = 0; v1 < number_of_vertices; v1++) {
		copy_tree.set_original_id(v1, tree.get_original_id(v1));
		for (int v2 = 0; v2 < number_of_vertices; v2++) {
			//std::cout<<"\nCopy Tree.  V1: " << v1 << "   V2: " << v2;
			if (tree.is_connected(v1,v2)) {
				//std::cout<<"\nTree is connected";
				for (int i = 0; i < tree.get_number_of_duplicate_edges(v1,v2); i++){ 
					copy_tree.addEdge(v1,v2);
				}
			}
		}
	}
	return copy_tree;
 }


/**
 * Given a vector of edges.  Checks if the given edge is in the vector
 */
bool AgentSolver::edge_in(std::vector<AgentSolver::Edge> edges, long id) {
	for (int i = 0; (unsigned) i < edges.size(); i++) {
		if (edges[i].getID() == id) {
			return true;
		}
	}
	return false;
}
/**
 * Is given the graph and two vertices.  It checks if there is a path between them
 */
bool AgentSolver::vertices_are_in_the_same_set(Adjacency_List al, int v1, int v2) {

	//List of all visited vertices
	std::vector<int> visited_vertices;
	
	//List of vertices in current search (i.e. they are all in the same set/partition)
	std::vector<int> vertice_set;
	
	int number_of_vertices = al.get_number_of_vertices();
	AgentSolver::Adjacency_List spanning_tree (number_of_vertices);
	
	//Used to keep track of where we are in graph
	std::stack<int> s;
	
	//Push first vertice onto stack
	s.push(0);
	visited_vertices.push_back(0);
	vertice_set.push_back(0);
	while(!all_vertices_visited(visited_vertices, al)) {
		while(!s.empty()) {
			int v1 = s.top();
			//Go through all vertices, 
			//check if it is connected to the current vertice 
			//and if it has not been visited before
			for(int v2 = 0; v2 < number_of_vertices; v2++) {
				if (al.is_connected(v1, v2) && std::find(visited_vertices.begin(), visited_vertices.end(), v2) == visited_vertices.end()) {
					spanning_tree.addEdge(v1,v2);
					spanning_tree.addEdge(v2,v1);
					s.push(v2);
					visited_vertices.push_back(v2);
					vertice_set.push_back(v2);
					v1 = s.top();
					v2 = 0;
				}
			}
			s.pop();
		}
		//If both vertices were visited in the same run
		if (std::find(vertice_set.begin(), vertice_set.end(), v1 ) != vertice_set.end() && std::find(vertice_set.begin(), vertice_set.end(), v2 ) != vertice_set.end()) {
			std::cout << "\nBoth Vertices visited in the same run";
			return true;
		}
		//Else if one or the other, but not both were found in the same run
		else if ((std::find(vertice_set.begin(), vertice_set.end(), v1 ) == vertice_set.end() && std::find(vertice_set.begin(), vertice_set.end(), v2 ) != vertice_set.end()) ||
				 (std::find(vertice_set.begin(), vertice_set.end(), v1 ) != vertice_set.end() && std::find(vertice_set.begin(), vertice_set.end(), v2 ) == vertice_set.end())){
			std::cout << "\nOne or the other visited in the same run.  Not both";
			return false;
		}		
		
		clear(s);
		
		//Clear the vertice set for the next search
		vertice_set.clear();
		
		
		//Pop the unvisited vertice onto the stack
		for(int v = 0; v < number_of_vertices; v++) {
			//If v is not in visited_vertices array, then it has not been visited
			//Push it onto stack and start depth first search again
			//If we find a vertice here, then we know the graph is not connected
			if (std::find(visited_vertices.begin(), visited_vertices.end(), v) == visited_vertices.end()) {
				s.push(v);
				visited_vertices.push_back(v);
				vertice_set.push_back(v);
				break;
			}
		}
	}
	return false;


}

/**
 * Is given the edges used so far of a graph, and the total graph
 * and determines if every edge of the total graph is in edges_used
 */
bool AgentSolver::not_all_edges_used(std::vector<AgentSolver::Edge> edges_used, Adjacency_List board_matrix) {
	if (edges_used.size() < (unsigned) board_matrix.get_number_of_edges()) {
		printf("\nUsed Edges Size: %d \nboard_matrix size %d", (signed)edges_used.size() ,board_matrix.get_number_of_edges());
	
		return true;
	}
	printf("\nUsed Edges Size: %d \nboard_matrix size %d", (signed)edges_used.size() ,board_matrix.get_number_of_edges());
	
	return false;
}


/**
 * Returns an adjacency list representing the board in its current state
 */ 
AgentSolver::Adjacency_List AgentSolver::getAdjacencyList() {	
	int size = rootboard.get_size();
	//Get the number of vertices on the bridg_it board 
	int numberOfVertices =  ((size - 1) / 2) * ((size - 1) / 2 - 1) + 2;
	//Get the size of the board

	//Create an adjacency list
	AgentSolver::Adjacency_List adjacency_list (numberOfVertices);

	//Get every connection
	//If 'Edge' contains opponent, then that edge is gone
	//If 'Edge' contains us, then that edge is reinforced
	Side edge;


	//printf("RootBoard.toplay() %d\n", rootboard.toplay().to_i());
	//Go through every piece on the board.  Look for vertices and check for connecting edges
	for (int xy = 0; xy < rootboard.vecsize(); xy++) {
		printf("XY: %d\n", xy);
		//Check if piece is a vertice of the person to play next
		if (rootboard.get(xy) == rootboard.toplay() && AgentSolver::xy_is_a_vertice(xy)) {
			printf("Got through\n");
			if (AgentSolver::xy_on_board(xy, xy + size * 2) && rootboard.get(xy + size) == Side::NONE) {
				printf("Piece right below is empty\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy + size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + size*2));
			}
			else if (AgentSolver::xy_on_board(xy, xy + size * 2) && rootboard.get(xy + size) == rootboard.toplay()) {
				printf("Piece right below is ours\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy + size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + size*2));
			}
		
			//Check directly above
			if (AgentSolver::xy_on_board(xy, xy - size * 2) && rootboard.get(xy - size) == Side::NONE) {
				printf("Piece right above is empty\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy - size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - size*2));
			}
			else if (AgentSolver::xy_on_board(xy, xy - size * 2) && rootboard.get(xy - size) == rootboard.toplay()) {
				printf("Piece right above is ours\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy - size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - size*2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - size*2));
			}
		
			//Check directly to the right
			if (AgentSolver::xy_on_board(xy, xy + 2) && rootboard.get(xy + 1) == Side::NONE) {
				printf("Piece to the right is empty\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy + 2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + 2));
			}
			else if (AgentSolver::xy_on_board(xy, xy + 2) && rootboard.get(xy + 1) == rootboard.toplay()) {
				printf("Piece to the right is ours\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy + 2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + 2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy + 2));
			}
			
			//Check directly to the left
			if (AgentSolver::xy_on_board(xy, xy - 2) && rootboard.get(xy - 1) == Side::NONE) {
				printf("Piece to the left is empty\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy - 2));
				adjacency_list.addEdge(xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - 2));
			}		
			else if (AgentSolver::xy_on_board(xy, xy - 2) && rootboard.get(xy - 1) == rootboard.toplay()) {
				printf("Piece to the left is ours\t");
				printf("Add edges: (%d, %d)\n", AgentSolver::xy_to_vertice(xy),  AgentSolver::xy_to_vertice(xy - 2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - 2));
				adjacency_list.addEdge(AgentSolver::xy_to_vertice(xy), AgentSolver::xy_to_vertice(xy - 2));
			}
		}
	}
	adjacency_list.addEdge(0, adjacency_list.get_number_of_vertices() - 1);
	adjacency_list.addEdge(adjacency_list.get_number_of_vertices() - 1, 0);
	adjacency_list.graph_to_s();
	return adjacency_list;
}

/**
 * Is given the xy value of the piece and returns which vertice # 
 * on the graph it represents
 */
int AgentSolver::xy_to_vertice(int xy) {
	if (rootboard.toplay() == Side::P1) {
		xy = AgentSolver::xy_from_whites_perspective(xy);
		}
	int size = rootboard.get_size();
	int vertice_number = 0;
	int i = 0;
	for (; i <= xy; i++) {
		Side vertice;
		if (rootboard.toplay() == Side::P1) {
			vertice = rootboard.get(xy_from_whites_perspective(i));
		}
		else {
			vertice = rootboard.get(i);
		}
		if ((i >= size) && (i < size * (size - 1)) && AgentSolver::xy_is_a_vertice(i) && (vertice == rootboard.toplay())) {
			vertice_number += 1;
		}
	}
	if (i >= (size * (size - 1))) {
		vertice_number += 1;
	}
	return vertice_number;
 }
 
 /**
  * Is given an edge, and returns an xy value representing the middle 
  * of that edge
  */
 int AgentSolver::edge_to_xy(int v1, int v2) {
	 int size = rootboard.get_size();
	 int final_xy = -1;
	 
	 
	 //Initialize XY of the vertices
	 int v1_xy = -1;
	 int v2_xy = -1;
	 std::cout<< "\n Edge V1: " << v1 << "  Edge V2: " << v2;
	for (int xy = rootboard.vecsize() - 1; xy >= 0 ; xy--) {
		 if (rootboard.toplay() == Side::P1) {
			 xy = xy_from_whites_perspective(xy);
		 }
		 std::cout << "\nxy_to_vertice(" << xy << ") = " << xy_to_vertice(xy); 
		 if (xy_to_vertice(xy) == v1) {
			 if (rootboard.toplay() == Side::P1) {
				v1_xy = xy_from_whites_perspective(xy);			 
			 }
			 else {
				v1_xy = xy;
			 }
		 }
		 else if (xy_to_vertice(xy) == v2) {
			 if (rootboard.toplay() == Side::P1) {
				v2_xy = xy_from_whites_perspective(xy);			 
			 }
			 else {
				v2_xy = xy;
			 }
		 }
		 if (rootboard.toplay() == Side::P1) {
			 xy = xy_from_whites_perspective(xy);
		 }		 
	 }
	 
	 
	 //Swap if v1_xy is on starting or finishing vertice
	 if (v1_xy == 0 || v1_xy >= (size*(size-1))) {
		 int temp = v1_xy;
		 v1_xy = v2_xy;
		 v2_xy = temp;
	 }
	 
	 std::cout<< "\nV1_XY = " << v1_xy << "\nV2_XY = " << v2_xy << "\n";
	 
		 
	 //4 Cases:
	 //If v2_xy is above v1_xy
	 if (v1_xy - v2_xy >= size) {
		 final_xy = v1_xy - size;
	 }
	 //If v2_xy is to the right
	 else if (v1_xy - v2_xy == -2) {
		 final_xy = v1_xy + 1;
	 }
	 //If v2_xy is below
	 else if (v1_xy - v2_xy <= -size) {
		 final_xy = v1_xy + size;
	 }
	 //If v2_xy is to the left
	 else if (v1_xy - v2_xy == 2) {
		 final_xy =  v1_xy - 1;
	 }
	
	if (rootboard.toplay() == Side::P1) {
		return xy_from_whites_perspective(final_xy);
	}
	else {
		return final_xy;
	}
	
}
 

bool AgentSolver::xy_is_a_vertice(int xy) {
	if (xy % 2 == 1) {
		return true;
	}
	else {
		return false;
	}
}

/**
 * Returns the total number of vertices on the board for one player
 */
int AgentSolver::get_number_of_vertices(Adjacency_List tree) {
	

	return tree.get_number_of_vertices();
}

/**
 * Returns the xy as if the board was mirrored on the diagonal line 
 * from top left to bottom right
 */
int AgentSolver::xy_from_whites_perspective(int xy) {
	int size = rootboard.get_size();
	int xy_from_whites_perspective = floor(xy / size) + (xy % size)*size;
	return xy_from_whites_perspective;
}

/**
 * Is given the xy of the original piece, and the connecting piece
 * Checks if the given connecting piece is on the board
 */
bool AgentSolver::xy_on_board(int xy, int xy2) {
	if (rootboard.get(xy) != rootboard.get(xy2)) {
		return false;
	}
	if (xy2 < 0 && xy2 >= rootboard.vecsize()) {
		return false;
	}
	return true;
}

void AgentSolver::set_board(const Board & board, bool clear){
	rootboard = board;
	Side a;
	
	AgentSolver::Adjacency_List board_matrix = AgentSolver::getAdjacencyList();
	board_matrix.graph_to_s();	
}


Move AgentSolver::return_move(const Node * node, Side toplay, int verbose) const {
	return root.bestmove;
}

long AgentSolver::get_id(){

    static long id = 0;
    id += 1;
    if (id > 1000000) {
		id = 0;
	}
    return id;
	}

long AgentSolver::get_random(long max) {
	srand(get_id());
	return (rand() % max);
}









































const float AgentSolver::min_rave = 0.1;

std::string AgentSolver::Node::to_s() const {
	return "AgentSolver::Node"
	       ", move " + move.to_s() +
	       ", exp " + exp.to_s() +
	       ", rave " + rave.to_s() +
	       ", know " + to_str(know) +
	       ", outcome " + to_str((int)outcome.to_i()) +
	       ", depth " + to_str((int)proofdepth) +
	       ", best " + bestmove.to_s() +
	       ", children " + to_str(children.num());
}

bool AgentSolver::Node::from_s(std::string s) {
	auto dict = parse_dict(s, ", ", " ");

	if(dict.size() == 9){
		move = Move(dict["move"]);
		exp = ExpPair(dict["exp"]);
		rave = ExpPair(dict["rave"]);
		know = from_str<int>(dict["know"]);
		outcome = Outcome(from_str<int>(dict["outcome"]));
		proofdepth = from_str<int>(dict["depth"]);
		bestmove = Move(dict["best"]);
		// ignore children
		return true;
	}
	return false;
}



AgentSolver::AgentSolver() : pool(this) {
	nodes = 0;
	runs = 0;
	gclimit = 5;

	profile     = false;
	ponder      = false;
//#ifdef SINGLE_THREAD ... make sure only 1 thread
	numthreads  = 1;
	pool.set_num_threads(numthreads);
	maxmem      = 1000*1024*1024;

	msrave      = -2;
	msexplore   = 0;

	explore     = 0;
	parentexplore = false;
	ravefactor  = 500;
	decrrave    = 0;
	knowledge   = true;
	userave     = 1;
	useexplore  = 1;
	fpurgency   = 1;
	rollouts    = 5;
	dynwiden    = 0;
	logdynwiden = (dynwiden ? std::log(dynwiden) : 0);

	shortrave   = false;
	keeptree    = true;
	minimax     = 2;
	visitexpand = 1;
	prunesymmetry = false;
	gcsolved    = 100000;

	localreply  = 5;
	locality    = 5;
	connect     = 20;
	size        = 0;
	bridge      = 100;
	dists       = 0;

	weightedrandom = false;
	rolloutpattern = true;
	lastgoodreply  = false;
	instantwin     = 0;

	for(int i = 0; i < 4096; i++)
		gammas[i] = 1;
}
AgentSolver::~AgentSolver(){
	pool.pause();
	pool.set_num_threads(0);

	root.dealloc(ctmem);
	ctmem.compact();
}

void AgentSolver::set_ponder(bool p){
	if(ponder != p){
		ponder = p;
		pool.pause();

		if(ponder)
			pool.resume();
	}
}


void AgentSolver::move(const Move & m){
	pool.pause();

	uword nodesbefore = nodes;

	if(keeptree && root.children.num() > 0){
		Node child;

		for(Node * i = root.children.begin(); i != root.children.end(); i++){
			if(i->move == m){
				child = *i;          //copy the child experience to temp
				child.swap_tree(*i); //move the child tree to temp
				break;
			}
		}

		nodes -= root.dealloc(ctmem);
		root = child;
		root.swap_tree(child);

		if(nodesbefore > 0)
			logerr("Nodes before: " + to_str(nodesbefore) + ", after: " + to_str(nodes) + ", saved " +  to_str(100.0*nodes/nodesbefore, 1) + "% of the tree\n");
	}else{
		nodes -= root.dealloc(ctmem);
		root = Node();
		root.move = m;
	}
	assert(nodes == root.size());

	rootboard.move(m);

	root.exp.addwins(visitexpand+1); //+1 to compensate for the virtual loss
	if(rootboard.won() < Outcome::DRAW)
		root.outcome = Outcome::UNKNOWN;

	if(ponder)
		pool.resume();
}

double AgentSolver::gamelen() const {
	DepthStats len;
	for(auto & t : pool)
		len += t->gamelen;
	return len.avg();
}

std::vector<Move> AgentSolver::get_pv() const {
	vecmove pv;

	const Node * n = & root;
	Side turn = rootboard.toplay();
	while(n && !n->children.empty()){
		Move m = return_move(n, turn);
		pv.push_back(m);
		n = find_child(n, m);
		turn = ~turn;
	}

	if(pv.size() == 0)
		pv.push_back(Move(M_RESIGN));

	return pv;
}

std::string AgentSolver::move_stats(vecmove moves) const {
	std::string s = "";
	const Node * node = & root;

	if(moves.size()){
		s += "path:\n";
		for(auto m : moves){
			if(node){
				node = find_child(node, m);
				s += node->to_s() + "\n";
			}
		}
	}

	if(node){
		s += "children:\n";
		for(auto & n : node->children)
			s += n.to_s() + "\n";
	}
	return s;
}


void AgentSolver::garbage_collect(Board & board, Node * node){
	Node * child = node->children.begin(),
		 * end = node->children.end();

	Side toplay = board.toplay();
	for( ; child != end; child++){
		if(child->children.num() == 0)
			continue;

		if(	(node->outcome >= Outcome::DRAW && child->exp.num() > gcsolved && (node->outcome != toplay || child->outcome == toplay || child->outcome == Outcome::DRAW)) || //parent is solved, only keep the proof tree, plus heavy draws
			(node->outcome <  Outcome::DRAW && child->exp.num() > (child->outcome >= Outcome::DRAW ? gcsolved : gclimit)) ){ // only keep heavy nodes, with different cutoffs for solved and unsolved
			board.set(child->move);
			garbage_collect(board, child);
			board.unset(child->move);
		}else{
			nodes -= child->dealloc(ctmem);
		}
	}
}

AgentSolver::Node * AgentSolver::find_child(const Node * node, const Move & move) const {
	for(auto & c : node->children)
		if(c.move == move)
			return &c;
	return NULL;
}

void AgentSolver::gen_sgf(SGFPrinter<Move> & sgf, unsigned int limit, const Node & node, Side side) const {
	for(auto & child : node.children){
		if(child.exp.num() >= limit && (side != node.outcome || child.outcome == node.outcome)){
			sgf.child_start();
			sgf.move(side, child.move);
			sgf.comment(child.to_s());
			gen_sgf(sgf, limit, child, ~side);
			sgf.child_end();
		}
	}
}

void AgentSolver::create_children_simple(const Board & board, Node * node){
	assert(node->children.empty());

	node->children.alloc(board.movesremain(), ctmem);

	Node * child = node->children.begin(),
		 * end   = node->children.end();
	Board::MoveIterator moveit = board.moveit(prunesymmetry);
	int nummoves = 0;
	for(; !moveit.done() && child != end; ++moveit, ++child){
		*child = Node(*moveit);
		nummoves++;
	}

	if(prunesymmetry)
		node->children.shrink(nummoves); //shrink the node to ignore the extra moves
	else //both end conditions should happen in parallel
		assert(moveit.done() && child == end);

	PLUS(nodes, node->children.num());
}

void AgentSolver::load_sgf(SGFParser<Move> & sgf, const Board & board, Node & node) {
	assert(sgf.has_children());
	create_children_simple(board, & node);

	while(sgf.next_child()){
		Move m = sgf.move();
		Node & child = *find_child(&node, m);
		child.from_s(sgf.comment());
		if(sgf.done_child()){
			continue;
		}else{
			// has children!
			Board b = board;
			b.move(m);
			load_sgf(sgf, b, child);
			assert(sgf.done_child());
		}
	}
}
































}; // namespace BridgIt
}; // namespace Morat
