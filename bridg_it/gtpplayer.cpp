

#include <fstream>

#include "../lib/fileio.h"

#include "gtp.h"

using namespace std;


GTPResponse GTP::gtp_move_stats(vecstr args){
	string s = "";

	Player::Node * node = &(player.root);

	for(unsigned int i = 0; i < args.size(); i++){
		Move m(args[i]);
		Player::Node * c = node->children.begin(),
		             * cend = node->children.end();
		for(; c != cend; c++){
			if(c->move == m){
				node = c;
				break;
			}
		}
	}

	Player::Node * child = node->children.begin(),
	             * childend = node->children.end();
	for( ; child != childend; child++){
		if(child->move == M_NONE)
			continue;

		s += child->move.to_s();
		s += "," + to_str((child->exp.num()  ? child->exp.avg() : 0.0), 4) + "," + to_str(child->exp.num());
		s += "," + to_str((child->rave.num() ? child->rave.avg() : 0.0), 4) + "," + to_str(child->rave.num());
		s += "," + to_str(child->know);
		if(child->outcome >= 0)
			s += "," + won_str(child->outcome);
		s += "\n";
	}
	return GTPResponse(true, s);
}

GTPResponse GTP::gtp_player_solve(vecstr args){
	double use_time = (args.size() >= 1 ?
			from_str<double>(args[0]) :
			time_control.get_time(hist.len(), hist->movesremain(), player.gamelen()));

	if(verbose)
		logerr("time remain: " + to_str(time_control.remain, 1) + ", time: " + to_str(use_time, 3) + ", sims: " + to_str(time_control.max_sims) + "\n");

	Player::Node * ret = player.genmove(use_time, time_control.max_sims, time_control.flexible);
	Move best = M_RESIGN;
	if(ret)
		best = ret->move;

	time_control.use(player.time_used);

	int toplay = player.rootboard.toplay();

	DepthStats gamelen, treelen;
	uint64_t runs = player.runs;
	double times[4] = {0,0,0,0};
	for(unsigned int i = 0; i < player.threads.size(); i++){
		gamelen += player.threads[i]->gamelen;
		treelen += player.threads[i]->treelen;

		for(int a = 0; a < 4; a++)
			times[a] += player.threads[i]->times[a];

		player.threads[i]->reset();
	}
	player.runs = 0;

	string stats = "Finished " + to_str(runs) + " runs in " + to_str(player.time_used*1000, 0) + " msec: " + to_str(runs/player.time_used, 0) + " Games/s\n";
	if(runs > 0){
		stats += "Game length: " + gamelen.to_s() + "\n";
		stats += "Tree depth:  " + treelen.to_s() + "\n";
		if(player.profile)
			stats += "Times:       " + to_str(times[0], 3) + ", " + to_str(times[1], 3) + ", " + to_str(times[2], 3) + ", " + to_str(times[3], 3) + "\n";
	}

	if(ret){
		stats += "Move Score:  " + to_str(ret->exp.avg()) + "\n";

		if(ret->outcome >= 0){
			stats += "Solved as a ";
			if(ret->outcome == toplay) stats += "win";
			else if(ret->outcome == 0) stats += "draw";
			else                       stats += "loss";
			stats += "\n";
		}
	}

	stats += "PV:          " + gtp_pv(vecstr()).response + "\n";

	if(verbose >= 3 && !player.root.children.empty())
		stats += "Exp-Rave:\n" + gtp_move_stats(vecstr()).response + "\n";

	if(verbose)
		logerr(stats);

	Solver s;
	if(ret){
		s.outcome = (ret->outcome >= 0 ? ret->outcome : -3);
		s.bestmove = ret->move;
		s.maxdepth = gamelen.maxdepth;
		s.nodes_seen = runs;
	}else{
		s.outcome = 3-toplay;
		s.bestmove = M_RESIGN;
		s.maxdepth = 0;
		s.nodes_seen = 0;
	}

	return GTPResponse(true, solve_str(s));
}


GTPResponse GTP::gtp_player_solved(vecstr args){
	string s = "";
	Player::Node * child = player.root.children.begin(),
	             * childend = player.root.children.end();
	int toplay = player.rootboard.toplay();
	int best = 0;
	for( ; child != childend; child++){
		if(child->move == M_NONE)
			continue;

		if(child->outcome == toplay)
			return GTPResponse(true, won_str(toplay));
		else if(child->outcome < 0)
			best = 2;
		else if(child->outcome == 0)
			best = 1;
	}
	if(best == 2) return GTPResponse(true, won_str(-3));
	if(best == 1) return GTPResponse(true, won_str(0));
	return GTPResponse(true, won_str(3 - toplay));
}

GTPResponse GTP::gtp_pv(vecstr args){
	string pvstr = "";
	vector<Move> pv = player.get_pv();
	for(unsigned int i = 0; i < pv.size(); i++)
		pvstr += pv[i].to_s() + " ";
	return GTPResponse(true, pvstr);
}

GTPResponse GTP::gtp_player_hgf(vecstr args){
	if(args.size() == 0)
		return GTPResponse(true, "player_hgf <filename> [sims limit]");

	FILE * fd = fopen(args[0].c_str(), "r");

	if(fd){
		fclose(fd);
		return GTPResponse(false, "File " + args[0] + " already exists");
	}

	fd = fopen(args[0].c_str(), "w");

	if(!fd)
		return GTPResponse(false, "Opening file " + args[0] + " for writing failed");

	unsigned int limit = 10000;
	if(args.size() > 1)
		limit = from_str<unsigned int>(args[1]);

	Board board = *hist;


	fprintf(fd, "(;FF[4]SZ[%i]\n", board.get_size());
	int p = 1;
	for(auto m : hist){
		fprintf(fd, ";%c[%s]", (p == 1 ? 'W' : 'B'), m.to_s().c_str());
		p = 3-p;
	}


	Player::Node * child = player.root.children.begin(),
	             * end = player.root.children.end();

	for( ; child != end; child++){
		if(child->exp.num() >= limit){
			board.set(child->move);
			player.gen_hgf(board, child, limit, 1, fd);
			board.unset(child->move);
		}
	}

	fprintf(fd, ")\n");

	fclose(fd);

	return true;
}

GTPResponse GTP::gtp_player_load_hgf(vecstr args){
	if(args.size() == 0)
		return GTPResponse(true, "player_load_hgf <filename>");

	FILE * fd = fopen(args[0].c_str(), "r");

	if(!fd)
		return GTPResponse(false, "Opening file " + args[0] + " for reading failed");

	int size;
	assert(fscanf(fd, "(;FF[4]SZ[%i]", & size) > 0);
	if(size != hist->get_size()){
		if(hist.len() == 0){
			hist = History(Board(size));
			set_board();
		}else{
			fclose(fd);
			return GTPResponse(false, "File has the wrong boardsize to match the existing game");
		}
	}

	eat_whitespace(fd);

	Board board(size);
	Player::Node * node = & player.root;
	vector<Player::Node *> prefix;

	char side, movestr[5];
	while(fscanf(fd, ";%c[%5[^]]]", &side, movestr) > 0){
		Move move(movestr);

		if(board.num_moves() >= (int)hist.len()){
			if(node->children.empty())
				player.create_children_simple(board, node);

			prefix.push_back(node);
			node = player.find_child(node, move);
		}else if(hist[board.num_moves()] != move){
			fclose(fd);
			return GTPResponse(false, "The current game is deeper than this file");
		}
		board.move(move);

		eat_whitespace(fd);
	}
	prefix.push_back(node);


	if(fpeek(fd) != ')'){
		if(node->children.empty())
			player.create_children_simple(board, node);

		while(fpeek(fd) != ')'){
			Player::Node child;
			player.load_hgf(board, & child, fd);

			Player::Node * i = player.find_child(node, child.move);
			*i = child;          //copy the child experience to the tree
			i->swap_tree(child); //move the child subtree to the tree

			assert(child.children.empty());

			eat_whitespace(fd);
		}
	}

	eat_whitespace(fd);
	assert(fgetc(fd) == ')');
	fclose(fd);

	while(!prefix.empty()){
		Player::Node * node = prefix.back();
		prefix.pop_back();

		Player::Node * child = node->children.begin(),
			         * end = node->children.end();

		int toplay = hist->toplay();
		if(prefix.size() % 2 == 1)
			toplay = 3 - toplay;

		Player::Node * backup = child;

		node->exp.clear();
		for( ; child != end; child++){
			node->exp += child->exp.invert();
			if(child->outcome == toplay || child->exp.num() > backup->exp.num())
				backup = child;
		}
		player.do_backup(node, backup, toplay);
	}

	return true;
}


GTPResponse GTP::gtp_genmove(vecstr args){
	if(player.rootboard.won() >= 0)
		return GTPResponse(true, "resign");

	double use_time = (args.size() >= 2 ?
			from_str<double>(args[1]) :
			time_control.get_time(hist.len(), hist->movesremain(), player.gamelen()));

	if(args.size() >= 2)
		use_time = from_str<double>(args[1]);

	if(verbose)
		logerr("time remain: " + to_str(time_control.remain, 1) + ", time: " + to_str(use_time, 3) + ", sims: " + to_str(time_control.max_sims) + "\n");

	uword nodesbefore = player.nodes;

	Player::Node * ret = player.genmove(use_time, time_control.max_sims, time_control.flexible);
	Move best = player.root.bestmove;

	time_control.use(player.time_used);

	int toplay = player.rootboard.toplay();

	DepthStats gamelen, treelen;
	uint64_t runs = player.runs;
	double times[4] = {0,0,0,0};
	for(unsigned int i = 0; i < player.threads.size(); i++){
		gamelen += player.threads[i]->gamelen;
		treelen += player.threads[i]->treelen;

		for(int a = 0; a < 4; a++)
			times[a] += player.threads[i]->times[a];

		player.threads[i]->reset();
	}
	player.runs = 0;

	string stats = "Finished " + to_str(runs) + " runs in " + to_str(player.time_used*1000, 0) + " msec: " + to_str(runs/player.time_used, 0) + " Games/s\n";
	if(runs > 0){
		stats += "Game length: " + gamelen.to_s() + "\n";
		stats += "Tree depth:  " + treelen.to_s() + "\n";
		if(player.profile)
			stats += "Times:       " + to_str(times[0], 3) + ", " + to_str(times[1], 3) + ", " + to_str(times[2], 3) + ", " + to_str(times[3], 3) + "\n";
	}

	if(ret)
		stats += "Move Score:  " + to_str(ret->exp.avg()) + "\n";

	if(player.root.outcome != -3){
		stats += "Solved as a ";
		if(player.root.outcome == 0)             stats += "draw";
		else if(player.root.outcome == toplay)   stats += "win";
		else if(player.root.outcome == 3-toplay) stats += "loss";
		else if(player.root.outcome == -toplay)  stats += "win or draw";
		else if(player.root.outcome == toplay-3) stats += "loss or draw";
		stats += "\n";
	}

	stats += "PV:          " + gtp_pv(vecstr()).response + "\n";

	if(verbose >= 3 && !player.root.children.empty())
		stats += "Exp-Rave:\n" + gtp_move_stats(vecstr()).response + "\n";

	string extended;
	if(genmoveextended){
		//move score
		if(ret) extended += " " + to_str(ret->exp.avg());
		else    extended += " 0";
		//outcome
		extended += " " + won_str(player.root.outcome);
		//work
		extended += " " + to_str(runs);
		//nodes
		extended += " " + to_str(player.nodes - nodesbefore);
	}

	move(best);

	if(verbose >= 2){
		stats += "history: ";
		for(auto m : hist)
			stats += m.to_s() + " ";
		stats += "\n";
		stats += hist->to_s(colorboard) + "\n";
	}

	if(verbose)
		logerr(stats);

	return GTPResponse(true, best.to_s() + extended);
}

GTPResponse GTP::gtp_player_params(vecstr args){
	if(args.size() == 0)
		return GTPResponse(true, string("\n") +
			"Set player parameters, eg: player_params -e 1 -f 0 -t 2 -o 1 -p 0\n" +
			"Processing:\n" +
#ifndef SINGLE_THREAD
			"  -t --threads     Number of MCTS threads                            [" + to_str(player.numthreads) + "]\n" +
#endif
			"  -o --ponder      Continue to ponder during the opponents time      [" + to_str(player.ponder) + "]\n" +
			"  -M --maxmem      Max memory in Mb to use for the tree              [" + to_str(player.maxmem/(1024*1024)) + "]\n" +
			"     --profile     Output the time used by each phase of MCTS        [" + to_str(player.profile) + "]\n" +
			"Final move selection:\n" +
			"  -E --msexplore   Lower bound constant in final move selection      [" + to_str(player.msexplore) + "]\n" +
			"  -F --msrave      Rave factor, 0 for pure exp, -1 # sims, -2 # wins [" + to_str(player.msrave) + "]\n" +
			"Tree traversal:\n" +
			"  -e --explore     Exploration rate for UCT                          [" + to_str(player.explore) + "]\n" +
			"  -A --parexplore  Multiply the explore rate by parents experience   [" + to_str(player.parentexplore) + "]\n" +
			"  -f --ravefactor  The rave factor: alpha = rf/(rf + visits)         [" + to_str(player.ravefactor) + "]\n" +
			"  -d --decrrave    Decrease the rave factor over time: rf += d*empty [" + to_str(player.decrrave) + "]\n" +
			"  -a --knowledge   Use knowledge: 0.01*know/sqrt(visits+1)           [" + to_str(player.knowledge) + "]\n" +
			"  -r --userave     Use rave with this probability [0-1]              [" + to_str(player.userave) + "]\n" +
			"  -X --useexplore  Use exploration with this probability [0-1]       [" + to_str(player.useexplore) + "]\n" +
			"  -u --fpurgency   Value to assign to an unplayed move               [" + to_str(player.fpurgency) + "]\n" +
			"  -O --rollouts    Number of rollouts to run per simulation          [" + to_str(player.rollouts) + "]\n" +
			"  -I --dynwiden    Dynamic widening, consider log_wid(exp) children  [" + to_str(player.dynwiden) + "]\n" +
			"Tree building:\n" +
			"  -s --shortrave   Only use moves from short rollouts for rave       [" + to_str(player.shortrave) + "]\n" +
			"  -k --keeptree    Keep the tree from the previous move              [" + to_str(player.keeptree) + "]\n" +
			"  -m --minimax     Backup the minimax proof in the UCT tree          [" + to_str(player.minimax) + "]\n" +
			"  -x --visitexpand Number of visits before expanding a node          [" + to_str(player.visitexpand) + "]\n" +
			"  -P --symmetry    Prune symmetric moves, good for proof, not play   [" + to_str(player.prunesymmetry) + "]\n" +
			"     --gcsolved    Garbage collect solved nodes with fewer sims than [" + to_str(player.gcsolved) + "]\n" +
			"Node initialization knowledge, Give a bonus:\n" +
			"  -l --localreply  based on the distance to the previous move        [" + to_str(player.localreply) + "]\n" +
			"  -y --locality    to stones near other stones of the same color     [" + to_str(player.locality) + "]\n" +
			"  -c --connect     to stones connected to edges                      [" + to_str(player.connect) + "]\n" +
			"  -S --size        based on the size of the group                    [" + to_str(player.size) + "]\n" +
			"  -b --bridge      to maintaining a 2-bridge after the op probes     [" + to_str(player.bridge) + "]\n" +
			"  -D --distance    to low minimum distance to win (<0 avoid VCs)     [" + to_str(player.dists) + "]\n" +
			"Rollout policy:\n" +
			"  -h --weightrand  Weight the moves according to computed gammas     [" + to_str(player.weightedrandom) + "]\n" +
			"  -p --pattern     Maintain the virtual connection pattern           [" + to_str(player.rolloutpattern) + "]\n" +
			"  -g --goodreply   Reuse the last good reply (1), remove losses (2)  [" + to_str(player.lastgoodreply) + "]\n" +
			"  -w --instantwin  Look for instant wins to this depth               [" + to_str(player.instantwin) + "]\n"
			);

	string errs;
	for(unsigned int i = 0; i < args.size(); i++) {
		string arg = args[i];

		if((arg == "-t" || arg == "--threads") && i+1 < args.size()){
			player.numthreads = from_str<int>(args[++i]);
			bool p = player.ponder;
			player.set_ponder(false); //stop the threads while resetting them
			player.reset_threads();
			player.set_ponder(p);
		}else if((arg == "-o" || arg == "--ponder") && i+1 < args.size()){
			player.set_ponder(from_str<bool>(args[++i]));
		}else if((arg == "--profile") && i+1 < args.size()){
			player.profile = from_str<bool>(args[++i]);
		}else if((arg == "-M" || arg == "--maxmem") && i+1 < args.size()){
			player.maxmem = from_str<uint64_t>(args[++i])*1024*1024;
		}else if((arg == "-E" || arg == "--msexplore") && i+1 < args.size()){
			player.msexplore = from_str<float>(args[++i]);
		}else if((arg == "-F" || arg == "--msrave") && i+1 < args.size()){
			player.msrave = from_str<float>(args[++i]);
		}else if((arg == "-e" || arg == "--explore") && i+1 < args.size()){
			player.explore = from_str<float>(args[++i]);
		}else if((arg == "-A" || arg == "--parexplore") && i+1 < args.size()){
			player.parentexplore = from_str<bool>(args[++i]);
		}else if((arg == "-f" || arg == "--ravefactor") && i+1 < args.size()){
			player.ravefactor = from_str<float>(args[++i]);
		}else if((arg == "-d" || arg == "--decrrave") && i+1 < args.size()){
			player.decrrave = from_str<float>(args[++i]);
		}else if((arg == "-a" || arg == "--knowledge") && i+1 < args.size()){
			player.knowledge = from_str<bool>(args[++i]);
		}else if((arg == "-s" || arg == "--shortrave") && i+1 < args.size()){
			player.shortrave = from_str<bool>(args[++i]);
		}else if((arg == "-k" || arg == "--keeptree") && i+1 < args.size()){
			player.keeptree = from_str<bool>(args[++i]);
		}else if((arg == "-m" || arg == "--minimax") && i+1 < args.size()){
			player.minimax = from_str<int>(args[++i]);
		}else if((arg == "-P" || arg == "--symmetry") && i+1 < args.size()){
			player.prunesymmetry = from_str<bool>(args[++i]);
		}else if((               arg == "--gcsolved") && i+1 < args.size()){
			player.gcsolved = from_str<uint>(args[++i]);
		}else if((arg == "-r" || arg == "--userave") && i+1 < args.size()){
			player.userave = from_str<float>(args[++i]);
		}else if((arg == "-X" || arg == "--useexplore") && i+1 < args.size()){
			player.useexplore = from_str<float>(args[++i]);
		}else if((arg == "-u" || arg == "--fpurgency") && i+1 < args.size()){
			player.fpurgency = from_str<float>(args[++i]);
		}else if((arg == "-O" || arg == "--rollouts") && i+1 < args.size()){
			player.rollouts = from_str<int>(args[++i]);
			if(player.gclimit < player.rollouts*5)
				player.gclimit = player.rollouts*5;
		}else if((arg == "-I" || arg == "--dynwiden") && i+1 < args.size()){
			player.dynwiden = from_str<float>(args[++i]);
			player.logdynwiden = std::log(player.dynwiden);
		}else if((arg == "-x" || arg == "--visitexpand") && i+1 < args.size()){
			player.visitexpand = from_str<uint>(args[++i]);
		}else if((arg == "-l" || arg == "--localreply") && i+1 < args.size()){
			player.localreply = from_str<int>(args[++i]);
		}else if((arg == "-y" || arg == "--locality") && i+1 < args.size()){
			player.locality = from_str<int>(args[++i]);
		}else if((arg == "-c" || arg == "--connect") && i+1 < args.size()){
			player.connect = from_str<int>(args[++i]);
		}else if((arg == "-S" || arg == "--size") && i+1 < args.size()){
			player.size = from_str<int>(args[++i]);
		}else if((arg == "-b" || arg == "--bridge") && i+1 < args.size()){
			player.bridge = from_str<int>(args[++i]);
		}else if((arg == "-D" || arg == "--distance") && i+1 < args.size()){
			player.dists = from_str<int>(args[++i]);
		}else if((arg == "-h" || arg == "--weightrand") && i+1 < args.size()){
			player.weightedrandom = from_str<int>(args[++i]);
		}else if((arg == "-p" || arg == "--pattern") && i+1 < args.size()){
			player.rolloutpattern = from_str<bool>(args[++i]);
		}else if((arg == "-g" || arg == "--goodreply") && i+1 < args.size()){
			player.lastgoodreply = from_str<int>(args[++i]);
		}else if((arg == "-w" || arg == "--instantwin") && i+1 < args.size()){
			player.instantwin = from_str<int>(args[++i]);
		}else{
			return GTPResponse(false, "Missing or unknown parameter");
		}
	}
	return GTPResponse(true, errs);
}

GTPResponse GTP::gtp_player_gammas(vecstr args){
	if(args.size() == 0)
		return GTPResponse(true, "Must pass the filename of a set of gammas");

	ifstream ifs(args[0].c_str());

	if(!ifs.good())
		return GTPResponse(false, "Failed to open file for reading");

	Board board = *hist;

	for(int i = 0; i < 4096; i++){
		int a;
		float f;
		ifs >> a >> f;

		if(i != a){
			ifs.close();
			return GTPResponse(false, "Line " + to_str(i) + " doesn't match the expected value");
		}

		int s = board.pattern_symmetry(i);
		if(s == i)
			player.gammas[i] = f;
		else
			player.gammas[i] = player.gammas[s];
	}

	ifs.close();
	return GTPResponse(true);
}
