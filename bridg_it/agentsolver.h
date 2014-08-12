
#pragma once


#include <cmath>
#include <cassert>
#include <list>
#include <stack>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "../lib/agentpool.h"
#include "../lib/compacttree.h"
#include "../lib/depthstats.h"
#include "../lib/exppair.h"
#include "../lib/log.h"
#include "../lib/thread.h"
#include "../lib/time.h"
#include "../lib/types.h"
#include "../lib/xorshift.h"



#include "agent.h"
#include "board.h"
#include "lbdist.h"
#include "move.h"
#include "movelist.h"
#include "policy_bridge.h"
#include "policy_instantwin.h"
#include "policy_lastgoodreply.h"
#include "policy_random.h"




namespace Morat {
namespace BridgIt {
	


class AgentSolver : public Agent{
public:

	struct Node {
	public:
		ExpPair rave;
		ExpPair exp;
		int16_t know;
		Outcome outcome;
		uint8_t proofdepth;
		Move    move;
		Move    bestmove; //if outcome is set, then bestmove is the way to get there
		CompactTree<Node>::Children children;
//		int padding;
		//seems to need padding to multiples of 8 bytes or it segfaults?
		//don't forget to update the copy constructor/operator

		Node() : know(0), outcome(Outcome::UNKNOWN), proofdepth(0), move(M_NONE) { }
		Node(const Move & m, Outcome o = Outcome::UNKNOWN) : know(0), outcome(o), proofdepth(0), move(m) { }
		Node(const Node & n) { *this = n; }
		Node & operator = (const Node & n){
			if(this != & n){ //don't copy to self
				//don't copy to a node that already has children
				assert(children.empty());

				rave = n.rave;
				exp  = n.exp;
				know = n.know;
				move = n.move;
				bestmove = n.bestmove;
				outcome = n.outcome;
				proofdepth = n.proofdepth;
				//children = n.children; ignore the children, they need to be swap_tree'd in
			}
			return *this;
		}

		void swap_tree(Node & n){
			children.swap(n.children);
		}

		std::string to_s() const ;
		bool from_s(std::string s);

		unsigned int size() const {
			unsigned int num = children.num();

			if(children.num())
				for(Node * i = children.begin(); i != children.end(); i++)
					num += i->size();

			return num;
		}

		~Node(){
			assert(children.empty());
		}

		unsigned int alloc(unsigned int num, CompactTree<Node> & ct){
			return children.alloc(num, ct);
		}
		unsigned int dealloc(CompactTree<Node> & ct){
			unsigned int num = 0;

			if(children.num())
				for(Node * i = children.begin(); i != children.end(); i++)
					num += i->dealloc(ct);
			num += children.dealloc(ct);

			return num;
		}

		//new way, more standard way of changing over from rave scores to real scores
		float value(float ravefactor, bool knowledge, float fpurgency){
			float val = fpurgency;
			float expnum = exp.num();
			float ravenum = rave.num();

			if(ravefactor <= min_rave){
				if(expnum > 0)
					val = exp.avg();
			}else if(ravenum > 0 || expnum > 0){
				float alpha = ravefactor/(ravefactor + expnum);
//				float alpha = sqrt(ravefactor/(ravefactor + 3.0f*expnum));
//				float alpha = ravenum/(expnum + ravenum + expnum*ravenum*ravefactor);

				val = 0;
				if(ravenum > 0) val += alpha*rave.avg();
				if(expnum  > 0) val += (1.0f-alpha)*exp.avg();
			}

			if(knowledge && know > 0){
				if(expnum <= 1)
					val += 0.01f * know;
				else if(expnum < 1000) //knowledge is only useful with little experience
					val += 0.01f * know / sqrt(expnum);
			}

			return val;
		}
	};
public:	

	struct Vertex {
		int id;
		int original_id;
		
		std::vector<int> attached;
		bool visited;
		bool exists;
		
		Vertex(int ID) {
			id = ID;
			original_id = ID;
			visited = false;
			exists = true;
		}
		Vertex(int ID, std::vector<int> Attached) {
			id = ID;
			original_id = ID;
			attached = Attached;
			visited = false;
			exists = true;
		}
		
		
		
		void add_attached(int v) {
			attached.push_back(v);
		}	
		int get_id() { 
			return id;
		}
		int get_original_id() {
			return original_id;
		} 
		
		void set_original_id(int id) {
			original_id = id;
		}
		
		int get_number_attached() {
			return attached.size();
		}
		std::vector<int> get_attached() {return attached;}
		bool get_visited() { return visited;}
		void set_visited(bool Visited) {visited = Visited;}
		void set_attached(std::vector<int> Attached) {attached = Attached;}
		void set_id(int ID) {id = ID;}
		void set_exists(bool Exists) {exists = Exists;}
		bool get_exists() {return exists;}
		bool is_attached(int v) {
			for (int i = 0; (unsigned)i < attached.size(); i++) {
				if(v == attached[i]) {
					return true;
				}
			}
			return false;
		}
		bool remove_attached(int v) {
			for (int i = 0; (unsigned)i < attached.size(); i++) {
				if (attached[i] == v) {
					attached.erase(attached.begin() + i);
					return true;
				}
			}
			return false;
		}
		
		void print() {
			for (int v = 0; (unsigned) v < attached.size(); v++) {
				std::cout << attached[v] << ", ";
			}
		}		
	};
	
	

	struct Edge {
		int v1;
		int v2; 
		long id;
		bool visited;
		
		Edge(int vert1, int vert2, long ID)  {
			v1 = vert1;
			v2 = vert2;
			id = ID;
			visited = false;
		}

		
		long getID() {
			return id;
		}
		
		int getV1() {
			return v1;
		}
		
		int getV2() {
			return v2;
		}
		
		void setV1(int v) {
			v1 = v;
		}
		
		void setV2(int v) {
			v2 = v;
		}
		
		bool get_visited() {
			return visited;
		}
		
	};
	

	class Adjacency_List {		
		private: 
			int size;
			std::vector<Vertex> vertices;
			std::vector<Edge> edges;
			
		public: 
			Adjacency_List() {}
			Adjacency_List(int size) {
				this->size = size;
				for(int i = 0; i < size; i++) {
					Vertex v = (i);
					vertices.push_back(v);
				}
			}
			
			
			void addEdge(Edge e) {
				int v1 = e.getV1();
				int v2 = e.getV2();
				
				addEdge(v1, v2, e.getID());
				addEdge(v2, v1, e.getID());				
			}
			
			void addEdge(int v1, int v2, int id=-1){
				if (v1 == v2){
					return;
				}
				//std::cout << "\nAdd Edge.  V1: " << v1 << "\tV2: " << v2<< "\n";
				for (int v = 0; (unsigned) v < vertices.size(); v++) {
					if (vertices[v].get_id() == v1) {
						v1 = v;
					}
				}
				vertices[v1].add_attached(v2);
				if (id == -1) {
					edges.push_back(Edge(v1, v2, edges.size()));
				}		
				else {
					edges.push_back(Edge(v1, v2, id));
				}					
			}
			
			std::vector<Edge> get_edges() {
				std::vector<Edge> culled_edges = edges; 
				for (int e = 0; (unsigned) e < culled_edges.size(); e++) {
					for (int i = e + 1; (unsigned)i < culled_edges.size(); i++) {
						if (culled_edges[e].getID() == culled_edges[i].getID()){
							culled_edges.erase(culled_edges.begin() + i);
						}			
					} 
				}
				return culled_edges;
			}
			
			int get_number_of_edges() {
				return edges.size()/2;
			}
			
			bool is_connected(Edge e1, Edge e2) {
				int e1v1 = e1.getV1();
				int e1v2 = e1.getV2();
				int e2v1 = e2.getV1();
				int e2v2 = e2.getV2();
				
				if (e1v1 == e2v2 || e1v1 == e2v1 ||
					e1v2 == e2v2 || e1v2 == e2v1) {
					return true;						
				}
				return false;
			}
			
			int get_number_of_edges_attached_to_vertice(int vertex) {
				for (int v = 0; (unsigned) v < vertices.size(); v++) {
					if (vertices[v].get_id() == vertex) {
						return vertices[v].get_number_attached();
					}
				}	
				std::cerr<<"Error!  Could not find specified vertice";			
				return -1;
			}
			
			int get_number_of_vertices() {
				return size;
			}
			
			int get_number_of_duplicate_edges(int v1, int v2) {
				int number = 0;
				for (int e = 0; (unsigned)e < edges.size(); e++) {
					if (edges[e].getV1() == v1 && edges[e].getV2() == v2) {
						number += 1;
					}
				}
				return number;
			}
			
			bool find_edge(int v1, int v2) {
				for(int e = 0; (unsigned)e < edges.size(); e++) {
					if (edges[e].getV1() == v1 && edges[e].getV2() == v2) {
						return true;
					}
				}
				return false;
			}
			
			bool find_edge(Edge edge) {
				for (int e = 0; (unsigned) e < edges.size(); e++) {
					if (edge.getID() == edges[e].getID()) {
						return true;
					}
				}
				return false;
			}
			
			
			bool delete_edge(int v1, int v2) {
				for (int v = 0; (unsigned) v < vertices.size(); v++) {
					if (vertices[v].get_id() == v1) {
						v1 = v;
						break;
					}
				}
				vertices[v1].remove_attached(v2);
				for(int e = 0; (unsigned)e < edges.size(); e++) {
					if (edges[e].getV1() == v1 && edges[e].getV2() == v2) {
						std::cout << "\nDeleting "<< v1 << " and " << v2;
						edges.erase(edges.begin() + e);
						return true;
					}
				}		
				return false;
			}
			
			
			bool delete_edge(Edge edge) {
				bool found_edge = false;
				int v1 = edge.getV1();
				int v2 = edge.getV2();
				for (int v = 0; (unsigned) v < vertices.size(); v++) {
					if (vertices[v].get_id() == v1) {
						v1 = v;
					}
				}
				vertices[v1].remove_attached(v2);
				vertices[v2].remove_attached(v1);
				for(int e = edges.size() - 1; e >= 0; e--) {
					if (edge.getID() == edges[e].getID()) {
						std::cout << "\nDeleting "<< v1 << " and " << v2;
						edges.erase(edges.begin() + e);
						found_edge = true;
					}
				}		
				return found_edge;
			}
			
			bool not_empty() {
				if (size != 0) {
					return true;
				}
				return false;
			}
			
			void set_original_id(int vertex, int original_id) {
				//Find the vertex specified and set its original id
				for (int v = 0; (unsigned)v < vertices.size(); v++) {
					if (vertices[v].get_id() == vertex) {
						vertices[v].set_original_id(original_id);
					}
				}
			}
			
			int get_original_id(int vertex) {
				for (int v = 0; (unsigned)v < vertices.size(); v++) {
					if (vertices[v].get_id() == vertex) {
						return vertices[v].get_original_id();
					}
				}
				return -1;
			}
			
			/**
			 * Deletes a vertice from the graph
			 */
			void delete_vertex(int vertex) {
				if (vertex >= size) {
					return;
				}	
				//Delete edges connected to this vertex
				for(int v = 0; (unsigned)v < vertices.size(); v++) {
					if (is_connected(v, vertex)) {
						delete_edge(v,vertex);
						delete_edge(vertex,v);
						v -= 1;
					}
				}
				//Delete all references of other vertices being attached to this vertex
				for(int v = 0; (unsigned)v < vertices.size(); v++) {
					if (is_connected(vertex, v)) {
						vertices[v].remove_attached(vertex);
					}
				}
				vertices.erase(vertices.begin() + vertex);
				size -= 1;
				

				
				//Reassign references for attached vertices
				for(int v = 0; (unsigned) v < vertices.size(); v++) {
					for (int r = vertex; (unsigned)r < vertices.size(); r++) { 
						if (is_connected(vertices[v].get_id(), vertices[r].get_id())) {
							std::cout<< "\nVertice: "<< vertices[v].get_id() << ", " << vertices[r].get_id() << "Are connected.";
							vertices[v].add_attached(vertices[r].get_id() - 1);
							vertices[v].remove_attached(vertices[r].get_id());
						}
					}
				}
				
				
				//Move all id's down one
				for (int v = vertex; (unsigned)v < vertices.size(); v++) {
					int id = vertices[v].get_id();
					vertices[v].set_id(id - 1);
				}
				
				for (int e = 0; (unsigned) e < edges.size(); e++) {
					if (edges[e].getV1() >= vertex) {
						edges[e].setV1(edges[e].getV1() - 1);
					}
					if (edges[e].getV2() >= vertex) {
						edges[e].setV2(edges[e].getV2() - 1);
					}					
				}
			}
			
				
			
						
			/**
			 * Returns true if the two are connected by a single edge
			 */
			bool is_connected(int v1, int v2) {
				for(int e = 0; (unsigned)e < edges.size(); e++) {
					if (edges[e].getV1() == v1 && edges[e].getV2() == v2) {
						return true;
					}
				}
				return false;
			}
			
			
			void graph_to_s(){        
				for (int v = 0; (unsigned)v < vertices.size(); v++) {
						std::cout << "\nVertex: " << vertices[v].get_id()<<  "(" << vertices[v].get_original_id()<< ")" << " : ";
						vertices[v].print();					
				}
			}

		
	};

	
	struct Vertice_List {
		std::vector<int> vertices;
		void addVertice(int vertice_number) {
			vertices.push_back(vertice_number);
		}
		std::vector<int> getVertices() {return vertices;}
		void print() {
			for (unsigned int i = 0; i < vertices.size(); i++) {
				std::cout << vertices[i] << ", ";
			}
		}
	};

	struct Partition {
		std::vector<Vertice_List> sets;
		Partition(int number_of_sets) : sets (number_of_sets) {}
		void addVertice(int vertice_number, int set_number) {
			sets[set_number].addVertice(vertice_number);
		}
		std::vector<int> getVertices(int set_number) {
			return sets[set_number].getVertices();
		}
		std::vector<Vertice_List> getSets() { return sets; }		
		void print() {
			std::cout << "\n****Partition*****\n";
			for(unsigned int i = 0; i < sets.size(); i++) {
				std::cout << "Set #" << i << ": ";
				sets[i].print();
				std::cout << "\n";
			}
		}
	};
	
	struct AdjacencyListNode{
		int destination;
		struct AdjacencyListNode *next;
	};


	struct AdjacencyList{
		struct AdjacencyListNode *head;
	};
	
	public:	
	Adjacency_List getAdjacencyList();
	int get_number_of_vertices(Adjacency_List tree);
	int xy_to_vertice(int xy);
	int xy_from_whites_perspective(int xy);
	bool xy_on_board(int xy, int xy2);
	bool xy_is_a_vertice(int xy);
	std::vector<Adjacency_List> find_edge_disjoint_trees(Adjacency_List board_matrix);
	std::vector<Adjacency_List> find_edge_disjoint_trees_old(Adjacency_List board_matrix);
	Adjacency_List get_spanning_tree(Adjacency_List board_matrix);
	bool is_connected(Adjacency_List tree);
	Adjacency_List subtract_trees(Adjacency_List tree1, Adjacency_List tree2);
	void swap_edges(Adjacency_List *tree1, Adjacency_List *tree2, std::vector<Edge> *used_edges, std::vector<Edge> *tree1_edges, std::vector<Edge> * tree2_edges);
	bool not_all_edges_used(std::vector<Edge> edges_used, Adjacency_List board_matrix);
	bool all_vertices_visited(std::vector<int> vertices_visited, Adjacency_List treey);
	void clear(std::stack<int> &s);
	bool vertices_are_in_the_same_set(Adjacency_List al, int v1, int v2);
	bool edge_in(std::vector<Edge> edges, long id);
	std::vector<Edge> get_all_edges(Adjacency_List tree);
	int get_id();
	Adjacency_List copyTree(Adjacency_List tree);
	long get_random(long max);
	Adjacency_List remove_problem_vertices(Adjacency_List tree);
	std::vector<Partition> push_to_all_indices(int i, std::vector<Partition> partitions);
	std::vector<Partition> concatenate(std::vector<Partition> p1, std::vector<Partition> p2);
	void get_best_move(std::vector<Adjacency_List> trees);
	int edge_to_xy(int v1, int v2);
	Move get_random_move();
	std::vector<Edge> get_cycle_edges(Adjacency_List tree, Edge e);
	void Augment( Adjacency_List * tree1, Adjacency_List * tree2, Adjacency_List * common_chords);

	  
	

	class AgentThread : public AgentThreadBase<AgentSolver> {
		mutable XORShift_float unitrand;
		LastGoodReply last_good_reply;
		RandomPolicy random_policy;
		ProtectBridge protect_bridge;
		InstantWin instant_wins;

		bool use_rave;    //whether to use rave for this simulation
		bool use_explore; //whether to use exploration for this simulation
		LBDists dists;    //holds the distances to the various non-ring wins as a heuristic for the minimum moves needed to win

		MoveList movelist;
		int stage; //which of the four MCTS stages is it on

	public:
		DepthStats treelen, gamelen;
		double times[4]; //time spent in each of the stages
		Time timestamps[4]; //timestamps for the beginning, before child creation, before rollout, after rollout

		AgentThread(AgentThreadPool<AgentSolver> * p, AgentSolver * a) : AgentThreadBase<AgentSolver>(p, a) { }


		void reset(){
			treelen.reset();
			gamelen.reset();

			use_rave = false;
			use_explore = false;

			for(int a = 0; a < 4; a++)
				times[a] = 0;
		}


	private:
		void iterate(); //handles each iteration
		void walk_tree(Board & board, Node * node, int depth);
		bool create_children(const Board & board, Node * node);
		void add_knowledge(const Board & board, Node * node, Node * child);
		Node * choose_move(const Node * node, Side toplay, int remain) const;
		void update_rave(const Node * node, Side toplay);
		bool test_bridge_probe(const Board & board, const Move & move, const Move & test) const;

		Outcome rollout(Board & board, Move move, int depth);
		Move rollout_choose_move(Board & board, const Move & prev);
		Move rollout_pattern(const Board & board, const Move & move);
	}; 

public:

	static const float min_rave;

	bool  ponder;     //think during opponents time?
	int   numthreads; //number of player threads to run
	u64   maxmem;     //maximum memory for the tree in bytes
	bool  profile;    //count how long is spent in each stage of MCTS
//final move selection
	float msrave;     //rave factor in final move selection, -1 means use number instead of value
	float msexplore;  //the UCT constant in final move selection
//tree traversal
	bool  parentexplore; // whether to multiple exploration by the parents winrate
	float explore;    //greater than one favours exploration, smaller than one favours exploitation
	float ravefactor; //big numbers favour rave scores, small ignore it
	float decrrave;   //decrease rave over time, add this value for each empty position on the board
	bool  knowledge;  //whether to include knowledge
	float userave;    //what probability to use rave
	float useexplore; //what probability to use UCT exploration
	float fpurgency;  //what value to return for a move that hasn't been played yet
	int   rollouts;   //number of rollouts to run after the tree traversal
	float dynwiden;   //dynamic widening, look at first log_dynwiden(experience) number of children, 0 to disable
	float logdynwiden; // = log(dynwiden), cached for performance
//tree building
	bool  shortrave;  //only update rave values on short rollouts
	bool  keeptree;   //reuse the tree from the previous move
	int   minimax;    //solve the minimax tree within the uct tree
	uint  visitexpand;//number of visits before expanding a node
	bool  prunesymmetry; //prune symmetric children from the move list, useful for proving but likely not for playing
	uint  gcsolved;   //garbage collect solved nodes or keep them in the tree, assuming they meet the required amount of work

//knowledge
	int   localreply; //boost for a local reply, ie a move near the previous move
	int   locality;   //boost for playing near previous stones
	int   connect;    //boost for having connections to edges and corners
	int   size;       //boost for large groups
	int   bridge;     //boost replying to a probe at a bridge
	int   dists;      //boost based on minimum number of stones needed to finish a non-ring win
//rollout
	int   weightedrandom; //use weighted random for move ordering based on gammas
	bool  rolloutpattern; //play the response to a virtual connection threat in rollouts
	int   lastgoodreply;  //use the last-good-reply rollout heuristic
	int   instantwin;     //how deep to look for instant wins in rollouts

	float gammas[4096]; //pattern weights for weighted random

	Board rootboard;
	Node  root;
	uword nodes;
	int   gclimit; //the minimum experience needed to not be garbage collected

	uint64_t runs, maxruns;

	CompactTree<Node> ctmem;

	AgentThreadPool<AgentSolver> pool;

	AgentSolver();
	~AgentSolver();

	void set_memlimit(uint64_t lim) { }; // in bytes
	void clear_mem() { };
	

	void set_ponder(bool p);
	void set_board(const Board & board, bool clear = true);

	void move(const Move & m);

	void search(double time, uint64_t maxruns, int verbose);
	Move return_move(int verbose) const { return return_move(& root, rootboard.toplay(), verbose); }

	double gamelen() const;
	vecmove get_pv() const;
	std::string move_stats(const vecmove moves) const;

	bool done() {
		//solved or finished runs
		return (rootboard.won() >= Outcome::DRAW || root.outcome >= Outcome::DRAW || (maxruns > 0 && runs >= maxruns));
	}

	bool need_gc() {
		//out of memory, start garbage collection
		return (ctmem.memalloced() >= maxmem);
	}

	void start_gc() {
		Time starttime;
		logerr("Starting player GC with limit " + to_str(gclimit) + " ... ");
		uint64_t nodesbefore = nodes;
		Board copy = rootboard;
		garbage_collect(copy, & root);
		Time gctime;
		ctmem.compact(1.0, 0.75);
		Time compacttime;
		logerr(to_str(100.0*nodes/nodesbefore, 1) + " % of tree remains - " +
			to_str((gctime - starttime)*1000, 0)  + " msec gc, " + to_str((compacttime - gctime)*1000, 0) + " msec compact\n");

		if(ctmem.meminuse() >= maxmem/2)
			gclimit = (int)(gclimit*1.3);
		else if(gclimit > rollouts*5)
			gclimit = (int)(gclimit*0.9); //slowly decay to a minimum of 5
	}

	void gen_sgf(SGFPrinter<Move> & sgf, int limit) const {
		if(limit < 0)
			limit = root.exp.num()/1000;
		gen_sgf(sgf, limit, root, rootboard.toplay());
	}

	void load_sgf(SGFParser<Move> & sgf) {
		load_sgf(sgf, rootboard, root);
	}

protected:

	void garbage_collect(Board & board, Node * node); //destroys the board, so pass in a copy
	bool do_backup(Node * node, Node * backup, Side toplay);
	Move return_move(const Node * node, Side toplay, int verbose = 0) const;

	Node * find_child(const Node * node, const Move & move) const ;
	void create_children_simple(const Board & board, Node * node);

	void gen_sgf(SGFPrinter<Move> & sgf, unsigned int limit, const Node & node, Side side) const ;
	void load_sgf(SGFParser<Move> & sgf, const Board & board, Node & node);


};






}; // namespace BridgIt
}; // namespace Morat
