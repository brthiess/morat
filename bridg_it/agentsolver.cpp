
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
namespace Hex {



/**
 * This function is called when genmove is called
 */
void AgentSolver::search(double time, uint64_t max_runs, int verbose){
	//Get the adjacency list representation of the board
	AgentSolver::Adjacency_List board_matrix = AgentSolver::getAdjacencyList();
	
	//Find out if we have won 
	find_winner(board_matrix);
	
	//Find 2 Edge Disjoint Spanning Trees
	find_edge_disjoint_trees(board_matrix);
}

/**
 * Is given an adjacency list of the current board state
 * 
 * Finds out if the person to play next has won the game
 * Uses the graph partitioning method
 */
Side AgentSolver::find_winner(Adjacency_List board_matrix) {
	
	//Creates partitions
	std::vector<AgentSolver::Partition> partitions = AgentSolver::get_partitions(board_matrix);
	
	//Simple Case:  Partition Graph into v sets where v = number of vertices
	int winning_number = 2 * (AgentSolver::get_number_of_vertices() - 1);
	
	//Check how many edges cross
	int number_of_edges = board_matrix.get_number_of_edges();
	
	//If the number of edges that cross is greater than the winning number
	//then we know there are at least 2 edge disjoint spanning trees
	//therefore we have won
	if (number_of_edges >= winning_number) {
		printf("Found win");
		return rootboard.toplay();
	}
	else if (number_of_edges >= winning_number - 1) {
		printf("Found win if playing first");
		return rootboard.toplay();
	}
	else {
		printf("Found a loss");
		return Side::NONE;
	}
	
}

/**
 * Is given a graph (represented by an adjacency list) 
 * And returns a vector of every partition
 */
 std::vector<AgentSolver::Partition> AgentSolver::get_partitions(Adjacency_List board_matrix) {
	 //Example vertices
	 int v1 = 1;
	 int v2 = 2;
	 int v3 = 3;
	 int v4 = 4;
	 
	 //Example Sets
	 int s0 = 0;
	 int s1 = 1;
	 int s2 = 2;


	 
	 //Create a partition with 3 sets
	 AgentSolver::Partition example_one (3);
	 
	 //Add vertice 4 to set 2, add vertice 3 to set 1, etc...
	 example_one.addVertice(v4, s2);
	 example_one.addVertice(v3, s1);
	 example_one.addVertice(v2, s1);
	 example_one.addVertice(v1, s0);
	 
	 //Print this partition to the terminal
	 example_one.print();
	 
	 //Create another example partition with 2 sets this time
	 AgentSolver::Partition example_two (2);
	 example_two.addVertice(v4, s0);
	 example_two.addVertice(v3, s1);
	 example_two.addVertice(v2, s1);
	 example_two.addVertice(v1, s0);	 
	 
	 //Print this partition to the terminal
	 example_two.print();
	 
	 
	 //Create a vector of example partitions
	 std::vector<AgentSolver::Partition> examples;
	 
	 //Add 'example_one' and 'example_two' to the examples vector
	 examples.push_back(example_one);
	 examples.push_back(example_two);
	 
	 //Return Examples
	 return examples;
 }

/**
 * Will return two edge disjoint trees, unless none are found in which case, 
 * it returns only one.
 */
std::vector<AgentSolver::Adjacency_List> AgentSolver::find_edge_disjoint_trees(Adjacency_List board_matrix) {
	
	//To keep track of the edges used for the swapping
	int number_of_vertices = AgentSolver::get_number_of_vertices();
	Adjacency_List edges_used (number_of_vertices);
	
	//Vector containing all edge disjoint trees found
	std::vector<AgentSolver::Adjacency_List> edge_disjoint_trees;
	
	//Find a spanning tree from the main graph
	Adjacency_List tree1 = get_spanning_tree(board_matrix);
	
	//Subtract the edges the edges from the main graph by tree1 and get a second tree
	Adjacency_List tree2 = subtract_trees(board_matrix, tree1);
	

	while(not_all_edges_used(edges_used, board_matrix)) {
		//Check if both trees are connected
		//If so then we have found 2 edge disjoint trees
		if (is_connected(tree2) && is_connected(tree1)) {
			edge_disjoint_trees.push_back(tree1);
			edge_disjoint_trees.push_back(tree2);
			break;
		}
		//Else find a connecting edge for tree2 (so that it is fully connected)
		//and swap it with tree1
		else {
			swap_edges(&tree1, &tree2, &edges_used);			
		}
	}
	
	
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
	
	int number_of_vertices = get_number_of_vertices();
	AgentSolver::Adjacency_List spanning_tree (number_of_vertices);
	
	//Used to keep track of where we are in graph
	std::stack<int> s;
	
	//Push first vertice onto stack
	s.push(0);
	visited_vertices.push_back(0);
	while(!all_vertices_visited(visited_vertices)) {
		while(!s.empty()) {
			int v1 = s.top();
			//Go through all vertices, 
			//check if it is connected to the current vertice 
			//and if it has not been visited before
			for(int v2 = 0; v2 < number_of_vertices; v2++) {
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
	spanning_tree.graph_to_s();
	return spanning_tree;
}

void AgentSolver::clear( std::stack<int> &s ) {
   std::stack<int> empty;
   std::swap( s, empty );
}


bool AgentSolver::all_vertices_visited(std::vector<int> visited_vertices) {
	int number_of_vertices = get_number_of_vertices();
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
AgentSolver::Adjacency_List AgentSolver::subtract_trees(Adjacency_List tree1, Adjacency_List tree2) {
	int number_of_vertices = get_number_of_vertices();
	AgentSolver::Adjacency_List subtracted_tree (number_of_vertices);
	//Go through all possible edge combinations and check to see if both graphs have them
	for (int v1 = 0; v1 < number_of_vertices; v1++) {
		for (int v2 = 0; v2 < number_of_vertices; v2++) {
			//If the main graph is connected and the second graph is not, then add it to the subtracted tree
			if (tree1.is_connected(v1,v2) && !tree2.is_connected(v1,v2)) {
				subtracted_tree.addEdge(v1,v2);
			}
		}	
	}
	subtracted_tree.graph_to_s();
	return subtracted_tree;
}

/**
 * Checks if the given tree connects every vertice
 * Returns true if it is connected
 */
bool AgentSolver::is_connected(Adjacency_List tree) {
	//List of all visited vertices
	std::vector<int> visited_vertices;
	
	int number_of_vertices = get_number_of_vertices();
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
		return false;
	}
	else {
		return true;
	}
}

/**
 * Is given two trees and a list of edges (all in pointer form) that have already been swapped with,
 * It finds an edge that will connect the unconnected tree
 * and that has not been swapped yet 
 * and swaps the two
 */
void AgentSolver::swap_edges(Adjacency_List *tree1, Adjacency_List *tree2, Adjacency_List *used_edges) {
	
}

/**
 * Is given the edges used so far of a graph, and the total graph
 * and determines if every edge of the total graph is in edges_used
 */
bool AgentSolver::not_all_edges_used(Adjacency_List edges_used, Adjacency_List board_matrix) {
	return false;
}


/**
 * Returns an adjacency list representing the board in its current state
 */ 
AgentSolver::Adjacency_List AgentSolver::getAdjacencyList() {	
	
	//Get the number of vertices on the bridg_it board 
	int numberOfVertices = AgentSolver::get_number_of_vertices();
	//Get the size of the board
	int size = rootboard.get_size();
	//Create an adjacency list
	AgentSolver::Adjacency_List adjacency_list (numberOfVertices);

	//Get every connection
	//If 'Edge' contains opponent, then that edge is gone
	//If 'Edge' contains us, then that edge is reinforced
	Side edge;


	printf("RootBoard.toplay() %d\n", rootboard.toplay().to_i());
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
 

bool AgentSolver::xy_is_a_vertice(int xy) {
	if (xy % 2 == 1 && AgentSolver::xy_from_whites_perspective(xy) % 2 == 1) {
		return true;
	}
	else {
		return false;
	}
}

/**
 * Returns the total number of vertices on the board for one player
 */
int AgentSolver::get_number_of_vertices() {
	int size = rootboard.get_size();

	int numberOfVertices = ((size - 1) / 2) * ((size - 1) / 2 - 1) + 2;

	return numberOfVertices;
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
	Move * m = new Move("a1");
	return *m;
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
































}; // namespace Hex
}; // namespace Morat
