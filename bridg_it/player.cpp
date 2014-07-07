
#include <cmath>
#include <string>

#include "../lib/alarm.h"
#include "../lib/fileio.h"
#include "../lib/string.h"
#include "../lib/time.h"

#include "board.h"
#include "player.h"

const float Player::min_rave = 0.1;

void Player::PlayerThread::run(){
	while(true){
		switch(player->threadstate){
		case Thread_Cancelled:  //threads should exit
			return;

		case Thread_Wait_Start: //threads are waiting to start
		case Thread_Wait_Start_Cancelled:
			player->runbarrier.wait();
			CAS(player->threadstate, Thread_Wait_Start, Thread_Running);
			CAS(player->threadstate, Thread_Wait_Start_Cancelled, Thread_Cancelled);
			break;

		case Thread_Wait_End:   //threads are waiting to end
			player->runbarrier.wait();
			CAS(player->threadstate, Thread_Wait_End, Thread_Wait_Start);
			break;

		case Thread_Running:    //threads are running
			if(player->rootboard.won() >= 0 || player->root.outcome >= 0 || (player->maxruns > 0 && player->runs >= player->maxruns)){ //solved or finished runs
				if(CAS(player->threadstate, Thread_Running, Thread_Wait_End) && player->root.outcome >= 0)
					logerr("Solved as " + to_str((int)player->root.outcome) + "\n");
				break;
			}
			if(player->ctmem.memalloced() >= player->maxmem){ //out of memory, start garbage collection
				CAS(player->threadstate, Thread_Running, Thread_GC);
				break;
			}

			INCR(player->runs);
			iterate();
			break;

		case Thread_GC:         //one thread is running garbage collection, the rest are waiting
		case Thread_GC_End:     //once done garbage collecting, go to wait_end instead of back to running
			if(player->gcbarrier.wait()){
				Time starttime;
				logerr("Starting player GC with limit " + to_str(player->gclimit) + " ... ");
				uint64_t nodesbefore = player->nodes;
				Board copy = player->rootboard;
				player->garbage_collect(copy, & player->root);
				Time gctime;
				player->ctmem.compact(1.0, 0.75);
				Time compacttime;
				logerr(to_str(100.0*player->nodes/nodesbefore, 1) + " % of tree remains - " +
					to_str((gctime - starttime)*1000, 0)  + " msec gc, " + to_str((compacttime - gctime)*1000, 0) + " msec compact\n");

				if(player->ctmem.meminuse() >= player->maxmem/2)
					player->gclimit = (int)(player->gclimit*1.3);
				else if(player->gclimit > player->rollouts*5)
					player->gclimit = (int)(player->gclimit*0.9); //slowly decay to a minimum of 5

				CAS(player->threadstate, Thread_GC,     Thread_Running);
				CAS(player->threadstate, Thread_GC_End, Thread_Wait_End);
			}
			player->gcbarrier.wait();
			break;
		}
	}
}

Player::Node * Player::genmove(double time, int max_runs, bool flexible){
	time_used = 0;
	int toplay = rootboard.toplay();

	if(rootboard.won() >= 0 || (time <= 0 && max_runs == 0))
		return NULL;

	Time starttime;

	stop_threads();

	if(runs)
		logerr("Pondered " + to_str(runs) + " runs\n");

	runs = 0;
	maxruns = max_runs;
	for(unsigned int i = 0; i < threads.size(); i++)
		threads[i]->reset();

	// if the move is forced and the time can be added to the clock, don't bother running at all
	if(!flexible || root.children.num() != 1){
		//let them run!
		start_threads();

		Alarm timer;
		if(time > 0)
			timer(time - (Time() - starttime), std::bind(&Player::timedout, this));

		//wait for the timer to stop them
		runbarrier.wait();
		CAS(threadstate, Thread_Wait_End, Thread_Wait_Start);
		assert(threadstate == Thread_Wait_Start);
	}

	if(ponder && root.outcome < 0)
		start_threads();

	time_used = Time() - starttime;

//return the best one
	return return_move(& root, toplay);
}



Player::Player() {
	nodes = 0;
	gclimit = 5;
	time_used = 0;

	profile     = false;
	ponder      = false;
//#ifdef SINGLE_THREAD ... make sure only 1 thread
	numthreads  = 1;
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

	weightedrandom = 0;
	rolloutpattern = true;
	lastgoodreply  = false;
	instantwin     = 0;

	for(int i = 0; i < 4096; i++)
		gammas[i] = 1;

	//no threads started until a board is set
	threadstate = Thread_Wait_Start;
}
Player::~Player(){
	stop_threads();

	numthreads = 0;
	reset_threads(); //shut down the theads properly

	root.dealloc(ctmem);
	ctmem.compact();
}
void Player::timedout() {
	CAS(threadstate, Thread_Running, Thread_Wait_End);
	CAS(threadstate, Thread_GC, Thread_GC_End);
}

string Player::statestring(){
	switch(threadstate){
	case Thread_Cancelled:  return "Thread_Wait_Cancelled";
	case Thread_Wait_Start: return "Thread_Wait_Start";
	case Thread_Wait_Start_Cancelled: return "Thread_Wait_Start_Cancelled";
	case Thread_Running:    return "Thread_Running";
	case Thread_GC:         return "Thread_GC";
	case Thread_GC_End:     return "Thread_GC_End";
	case Thread_Wait_End:   return "Thread_Wait_End";
	}
	return "Thread_State_Unknown!!!";
}

void Player::stop_threads(){
	if(threadstate != Thread_Wait_Start){
		timedout();
		runbarrier.wait();
		CAS(threadstate, Thread_Wait_End, Thread_Wait_Start);
		assert(threadstate == Thread_Wait_Start);
	}
}

void Player::start_threads(){
	assert(threadstate == Thread_Wait_Start);
	runbarrier.wait();
	CAS(threadstate, Thread_Wait_Start, Thread_Running);
}

void Player::reset_threads(){ //start and end with threadstate = Thread_Wait_Start
	assert(threadstate == Thread_Wait_Start);

//wait for them to all get to the barrier
	assert(CAS(threadstate, Thread_Wait_Start, Thread_Wait_Start_Cancelled));
	runbarrier.wait();

//make sure they exited cleanly
	for(unsigned int i = 0; i < threads.size(); i++){
		threads[i]->join();
		delete threads[i];
	}

	threads.clear();

	threadstate = Thread_Wait_Start;

	runbarrier.reset(numthreads + 1);
	gcbarrier.reset(numthreads);

//start new threads
	for(int i = 0; i < numthreads; i++)
		threads.push_back(new PlayerUCT(this));
}

void Player::set_ponder(bool p){
	if(ponder != p){
		ponder = p;
		stop_threads();

		if(ponder)
			start_threads();
	}
}

void Player::set_board(const Board & board){
	stop_threads();

	nodes -= root.dealloc(ctmem);
	root = Node();
	root.exp.addwins(visitexpand+1);

	rootboard = board;

	reset_threads(); //needed since the threads aren't started before a board it set

	if(ponder)
		start_threads();
}
void Player::move(const Move & m){
	stop_threads();

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
	if(rootboard.won() < 0)
		root.outcome = -3;

	if(ponder)
		start_threads();
}

double Player::gamelen(){
	DepthStats len;
	for(unsigned int i = 0; i < threads.size(); i++)
		len += threads[i]->gamelen;
	return len.avg();
}

vector<Move> Player::get_pv(){
	vector<Move> pv;

	Node * r, * n = & root;
	char turn = rootboard.toplay();
	while(!n->children.empty()){
		r = return_move(n, turn);
		if(!r) break;
		pv.push_back(r->move);
		turn = 3 - turn;
		n = r;
	}

	if(pv.size() == 0)
		pv.push_back(Move(M_RESIGN));

	return pv;
}

Player::Node * Player::return_move(Node * node, int toplay) const {
	double val, maxval = -1000000000000.0; //1 trillion

	Node * ret = NULL,
		 * child = node->children.begin(),
		 * end = node->children.end();

	for( ; child != end; child++){
		if(child->outcome >= 0){
			if(child->outcome == toplay) val =  800000000000.0 - child->exp.num(); //shortest win
			else if(child->outcome == 0) val = -400000000000.0 + child->exp.num(); //longest tie
			else                         val = -800000000000.0 + child->exp.num(); //longest loss
		}else{ //not proven
			if(msrave == -1) //num simulations
				val = child->exp.num();
			else if(msrave == -2) //num wins
				val = child->exp.sum();
			else
				val = child->value(msrave, 0, 0) - msexplore*sqrt(log(node->exp.num())/(child->exp.num() + 1));
		}

		if(maxval < val){
			maxval = val;
			ret = child;
		}
	}

//set bestmove, but don't touch outcome, if it's solved that will already be set, otherwise it shouldn't be set
	if(ret){
		node->bestmove = ret->move;
	}else if(node->bestmove == M_UNKNOWN){
		// TODO: Is this needed?
//		SolverAB solver;
//		solver.set_board(rootboard);
//		solver.solve(0.1);
//		node->bestmove = solver.bestmove;
	}

	assert(node->bestmove != M_UNKNOWN);

	return ret;
}

void Player::garbage_collect(Board & board, Node * node){
	Node * child = node->children.begin(),
		 * end = node->children.end();

	int toplay = board.toplay();
	for( ; child != end; child++){
		if(child->children.num() == 0)
			continue;

		if(	(node->outcome >= 0 && child->exp.num() > gcsolved && (node->outcome != toplay || child->outcome == toplay || child->outcome == 0)) || //parent is solved, only keep the proof tree, plus heavy draws
			(node->outcome <  0 && child->exp.num() > (child->outcome >= 0 ? gcsolved : gclimit)) ){ // only keep heavy nodes, with different cutoffs for solved and unsolved
			board.set(child->move);
			garbage_collect(board, child);
			board.unset(child->move);
		}else{
			nodes -= child->dealloc(ctmem);
		}
	}
}

Player::Node * Player::find_child(Node * node, const Move & move){
	for(Node * i = node->children.begin(); i != node->children.end(); i++)
		if(i->move == move)
			return i;

	return NULL;
}

void Player::gen_hgf(Board & board, Node * node, unsigned int limit, unsigned int depth, FILE * fd){
	string s = string("\n") + string(depth, ' ') + "(;" + (board.toplay() == 2 ? "W" : "B") + "[" + node->move.to_s() + "]" +
	       "C[mcts, sims:" + to_str(node->exp.num()) + ", avg:" + to_str(node->exp.avg(), 4) + ", outcome:" + to_str((int)(node->outcome)) + ", best:" + node->bestmove.to_s() + "]";
	fprintf(fd, "%s", s.c_str());

	Node * child = node->children.begin(),
		 * end = node->children.end();

	int toplay = board.toplay();

	bool children = false;
	for( ; child != end; child++){
		if(child->exp.num() >= limit && (toplay != node->outcome || child->outcome == node->outcome) ){
			board.set(child->move);
			gen_hgf(board, child, limit, depth+1, fd);
			board.unset(child->move);
			children = true;
		}
	}

	if(children)
		fprintf(fd, "\n%s", string(depth, ' ').c_str());
	fprintf(fd, ")");
}

void Player::create_children_simple(const Board & board, Node * node){
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

//reads the format from gen_hgf.
void Player::load_hgf(Board board, Node * node, FILE * fd){
	char c, buf[101];

	eat_whitespace(fd);

	assert(fscanf(fd, "(;%c[%100[^]]]", &c, buf) > 0);

	assert(board.toplay() == (c == 'W' ? 1 : 2));
	node->move = Move(buf);
	board.move(node->move);

	assert(fscanf(fd, "C[%100[^]]]", buf) > 0);

	vecstr entry, parts = explode(string(buf), ", ");
	assert(parts[0] == "mcts");

	entry = explode(parts[1], ":");
	assert(entry[0] == "sims");
	uword sims = from_str<uword>(entry[1]);

	entry = explode(parts[2], ":");
	assert(entry[0] == "avg");
	double avg = from_str<double>(entry[1]);

	uword wins = sims*avg;
	node->exp.addwins(wins);
	node->exp.addlosses(sims - wins);

	entry = explode(parts[3], ":");
	assert(entry[0] == "outcome");
	node->outcome = from_str<int>(entry[1]);

	entry = explode(parts[4], ":");
	assert(entry[0] == "best");
	node->bestmove = Move(entry[1]);


	eat_whitespace(fd);

	if(fpeek(fd) != ')'){
		create_children_simple(board, node);

		while(fpeek(fd) != ')'){
			Node child;
			load_hgf(board, & child, fd);

			Node * i = find_child(node, child.move);
			*i = child;          //copy the child experience to the tree
			i->swap_tree(child); //move the child subtree to the tree

			assert(child.children.empty());

			eat_whitespace(fd);
		}
	}

	eat_char(fd, ')');

	return;
}
