
#pragma once

//Interface for the various agents: players and solvers

#include "../lib/types.h"

#include "board.h"

class Agent {
public:
	Agent() { }
	virtual ~Agent() { }

	virtual void search(double time, uint64_t maxruns, int verbose) = 0;
	virtual Move return_move(int verbose) const = 0;
	virtual void set_board(const Board & board, bool clear = true) = 0;
	virtual void move(const Move & m) = 0;
	virtual void set_memlimit(uint64_t lim) = 0; // in bytes
	virtual void clear_mem() = 0;

	virtual vector<Move> get_pv() const = 0;
	        string move_stats() const { return move_stats(vector<Move>()); }
	virtual string move_stats(const vector<Move> moves) const = 0;
	virtual double gamelen() const = 0;

	virtual void timedout(){ timeout = true; }

protected:
	volatile bool timeout;
	Board rootboard;

	static int solve1ply(const Board & board, unsigned int & nodes) {
		int outcome = -3;
		int turn = board.toplay();
		for(Board::MoveIterator move = board.moveit(true); !move.done(); ++move){
			++nodes;
			int won = board.test_win(*move, turn);

			if(won == turn)
				return won;
			if(won == 0)
				outcome = 0;
		}
		return outcome;
	}

	static int solve2ply(const Board & board, unsigned int & nodes) {
		int losses = 0;
		int outcome = -3;
		int turn = board.toplay(), opponent = 3 - turn;
		if (board.won() > 0) {
			//printf("Game is already over.  Don't search here!");
			return board.won();
		}
		for(Board::MoveIterator move = board.moveit(true); !move.done(); ++move){
			++nodes;
			int won = board.test_win(*move, turn);
			
			/*printf("\n\nMove x: %d y: %d\n", move->x, move->y);
			board.to_s(true);
			board.print(true);*/

			if(won == turn) {
				//printf("Our Win\n");
				return won;
			}
			if(won == 0) {
				//printf("Outcome is 0?\n");
				outcome = 0;
			}

			if(board.test_win(*move, opponent) > 0) {
				//printf("Opponent Wins playing here\n");
				losses++;
			}
			//printf("Outcome: %d\n", outcome);
		}
		if(losses >= 2) {
			//printf("More than 2 losses");
			return opponent;
		}
		//printf("Return Outcome %d\n", outcome);
		return outcome;
	}

};
