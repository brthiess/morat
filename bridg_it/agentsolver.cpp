
#include <cmath>
#include <string>

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
 * 				1		First Row
 * 			  / | \
 * 			 2--3--4	Second Row
 * 			 |  |  |
 * 			 5--6--7	Third Row
 *  		  \ | /
 * 				8		Last Row 
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
}

/**
 * Is given an adjacency list of the current board state
 * 
 * Finds out if the person to play next has won the game
 * Uses the graph partitioning method
 */
Side AgentSolver::find_winner(Adjacency_List board_matrix) {
	
	//Simple Case:  Partition Graph into v sets where v = number of vertices
	int winning_number = 2 * (AgentSolver::get_number_of_vertices() - 1);
	
	//Check how many edges cross
	int number_of_edges = board_matrix.get_number_of_edges();
	
	//If the number of edges that cross is greater than the winning number
	//then we know there are at least 2 edge disjoint spanning trees
	//therefore we have won
	if (number_of_edges >= winning_number) {
		return rootboard.toplay();
	}
	else {
		return Side::NONE;
	}
	
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
		if ((i >= size) && (i < size * (size - 1)) && AgentSolver::xy_is_a_vertice(xy) && (vertice == rootboard.toplay())) {
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
