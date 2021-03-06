
#pragma once

#include "../lib/gtpcommon.h"
#include "../lib/string.h"

#include "agent.h"
//#include "agentab.h"
#include "agentmcts.h"
#include "agentpns.h"
#include "board.h"
#include "history.h"
#include "move.h"


namespace Morat {
namespace Rex {

class GTP : public GTPCommon {
	History hist;

public:
	int verbose;
	bool genmoveextended;
	bool colorboard;

	int mem_allowed;

	Agent * agent;

	GTP(FILE * i = stdin, FILE * o = stdout) : GTPCommon(i, o), hist(Board(Board::default_size)) {
		verbose = 1;
		genmoveextended = false;
		colorboard = true;

		mem_allowed = 1000;

		agent = new AgentMCTS();

		set_board();

		newcallback("name",            std::bind(&GTP::gtp_name,          this, _1), "Name of the program");
		newcallback("version",         std::bind(&GTP::gtp_version,       this, _1), "Version of the program");
		newcallback("verbose",         std::bind(&GTP::gtp_verbose,       this, _1), "Set verbosity, 0 for quiet, 1 for normal, 2+ for more output");
		newcallback("extended",        std::bind(&GTP::gtp_extended,      this, _1), "Output extra stats from genmove in the response");
		newcallback("debug",           std::bind(&GTP::gtp_debug,         this, _1), "Enable debug mode");
		newcallback("colorboard",      std::bind(&GTP::gtp_colorboard,    this, _1), "Turn on or off the colored board");
		newcallback("showboard",       std::bind(&GTP::gtp_print,         this, _1), "Show the board");
		newcallback("print",           std::bind(&GTP::gtp_print,         this, _1), "Alias for showboard");
		newcallback("dists",           std::bind(&GTP::gtp_dists,         this, _1), "Similar to print, but shows minimum win distances");
		newcallback("zobrist",         std::bind(&GTP::gtp_zobrist,       this, _1), "Output the zobrist hash for the current move");
		newcallback("clear_board",     std::bind(&GTP::gtp_clearboard,    this, _1), "Clear the board, but keep the size");
		newcallback("clear",           std::bind(&GTP::gtp_clearboard,    this, _1), "Alias for clear_board");
		newcallback("boardsize",       std::bind(&GTP::gtp_boardsize,     this, _1), "Clear the board, set the board size");
		newcallback("size",            std::bind(&GTP::gtp_boardsize,     this, _1), "Alias for board_size");
		newcallback("play",            std::bind(&GTP::gtp_play,          this, _1), "Place a stone: play <color> <location>");
		newcallback("white",           std::bind(&GTP::gtp_playwhite,     this, _1), "Place a white stone: white <location>");
		newcallback("black",           std::bind(&GTP::gtp_playblack,     this, _1), "Place a black stone: black <location>");
		newcallback("undo",            std::bind(&GTP::gtp_undo,          this, _1), "Undo one or more moves: undo [amount to undo]");
		newcallback("time",            std::bind(&GTP::gtp_time,          this, _1), "Set the time limits and the algorithm for per game time");
		newcallback("genmove",         std::bind(&GTP::gtp_genmove,       this, _1), "Generate a move: genmove [color] [time]");
		newcallback("solve",           std::bind(&GTP::gtp_solve,         this, _1), "Try to solve this position");

//		newcallback("ab",              std::bind(&GTP::gtp_ab,            this, _1), "Switch to use the Alpha/Beta agent to play/solve");
		newcallback("mcts",            std::bind(&GTP::gtp_mcts,          this, _1), "Switch to use the Monte Carlo Tree Search agent to play/solve");
		newcallback("pns",             std::bind(&GTP::gtp_pns,           this, _1), "Switch to use the Proof Number Search agent to play/solve");

		newcallback("all_legal",       std::bind(&GTP::gtp_all_legal,     this, _1), "List all legal moves");
		newcallback("history",         std::bind(&GTP::gtp_history,       this, _1), "List of played moves");
		newcallback("playgame",        std::bind(&GTP::gtp_playgame,      this, _1), "Play a list of moves");
		newcallback("winner",          std::bind(&GTP::gtp_winner,        this, _1), "Check the winner of the game");
		newcallback("patterns",        std::bind(&GTP::gtp_patterns,      this, _1), "List all legal moves plus their local pattern");

		newcallback("pv",              std::bind(&GTP::gtp_pv,            this, _1), "Output the principle variation for the player tree as it stands now");
		newcallback("move_stats",      std::bind(&GTP::gtp_move_stats,    this, _1), "Output the move stats for the player tree as it stands now");

		newcallback("params",          std::bind(&GTP::gtp_params,        this, _1), "Set the options for the player, no args gives options");

		newcallback("save_sgf",        std::bind(&GTP::gtp_save_sgf,      this, _1), "Output an sgf of the current tree");
		newcallback("load_sgf",        std::bind(&GTP::gtp_load_sgf,      this, _1), "Load an sgf generated by save_sgf");
//		newcallback("player_gammas",   std::bind(&GTP::gtp_player_gammas, this, _1), "Load the gammas for weighted random from a file");
	}

	void set_board(bool clear = true){
		agent->set_board(*hist);
	}

	void move(const Move & m){
		hist.move(m);
		agent->move(m);
	}

	GTPResponse gtp_print(vecstr args);
	GTPResponse gtp_zobrist(vecstr args);
	GTPResponse gtp_boardsize(vecstr args);
	GTPResponse gtp_clearboard(vecstr args);
	GTPResponse gtp_undo(vecstr args);
	GTPResponse gtp_all_legal(vecstr args);
	GTPResponse gtp_history(vecstr args);
	GTPResponse gtp_patterns(vecstr args);
	GTPResponse play(const std::string & pos, Side toplay);
	GTPResponse gtp_playgame(vecstr args);
	GTPResponse gtp_play(vecstr args);
	GTPResponse gtp_playwhite(vecstr args);
	GTPResponse gtp_playblack(vecstr args);
	GTPResponse gtp_winner(vecstr args);
	GTPResponse gtp_name(vecstr args);
	GTPResponse gtp_version(vecstr args);
	GTPResponse gtp_verbose(vecstr args);
	GTPResponse gtp_extended(vecstr args);
	GTPResponse gtp_colorboard(vecstr args);
	GTPResponse gtp_debug(vecstr args);
	GTPResponse gtp_dists(vecstr args);


	GTPResponse gtp_move_stats(vecstr args);
	GTPResponse gtp_pv(vecstr args);
	GTPResponse gtp_genmove(vecstr args);
	GTPResponse gtp_solve(vecstr args);

	GTPResponse gtp_params(vecstr args);

//	GTPResponse gtp_ab(vecstr args);
//	GTPResponse gtp_ab_params(vecstr args);
	GTPResponse gtp_mcts(vecstr args);
	GTPResponse gtp_mcts_params(vecstr args);
	GTPResponse gtp_pns(vecstr args);
	GTPResponse gtp_pns_params(vecstr args);

//	GTPResponse gtp_player_gammas(vecstr args);
	GTPResponse gtp_save_sgf(vecstr args);
	GTPResponse gtp_load_sgf(vecstr args);

	std::string solve_str(int outcome) const;
};

}; // namespace Rex
}; // namespace Morat
