
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
		
		//int losses = 0;
		int outcome = -3;
		int numberOfSelfToxicCells = 0;
		int numberOfOpponentToxicCells = 0;
		int turn = board.toplay(), opponent = 3 - turn;
		if (board.won() > 0) {
			//printf("Game is already over.  Don't search here!");
			return board.won();
		}
		for(Board::MoveIterator move = board.moveit(true); !move.done(); ++move){
			++nodes;
			int won = board.test_win(*move, turn);
			int opponent_loss = board.test_win(*move, opponent);
			
			/*printf("\n\nMove x: %d y: %d\n", move->x, move->y);
			board.to_s(true);
			board.print(true);*/
			//If opponent loses playing this move
			if(opponent_loss == turn) {
				//printf("Opponent Loss\n");
				numberOfOpponentToxicCells++;
				//proven loss for opponent
				if (numberOfOpponentToxicCells >= 2 || (numberOfOpponentToxicCells >= 1 && board.movesremain() % 2 == 0))
					return turn;
			}
			//If we lose playing this move
			else if(won == opponent) {
				//printf("Our loss\n");
				numberOfSelfToxicCells++;	
				//Proven Loss for us
				/*if(numberOfSelfToxicCells >= 2 || (board.movesremain() % 2 == 1 && numberOfSelfToxicCells >= 1))
					return opponent;*/
			}
			
			//Proven Win
			if(board.movesremain() == 2 && won < 0) {
				//printf("Proven win\n");
				return turn;
			}
		}
		//Proven Loss
		if(numberOfSelfToxicCells == board.movesremain()) {
			/*board.to_s(true);
			board.print(true);*/
			//printf("Board.movesremain() = %d\t Proven Loss\n", board.movesremain());
			return opponent;
		}
		/*
		if(losses >= 2)
			return opponent;*/
			//printf("Outcome: %d\n", outcome); 
		return outcome;
	}

};
