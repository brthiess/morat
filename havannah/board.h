
#pragma once

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

#include "../lib/bitcount.h"
#include "../lib/hashset.h"
#include "../lib/outcome.h"
#include "../lib/string.h"
#include "../lib/types.h"
#include "../lib/zobrist.h"

#include "move.h"

namespace Morat {
namespace Havannah {

/*
 * the board is represented as a flattened 2d array of the form:
 *   1 2 3
 * A 0 1 2    0 1       0 1
 * B 3 4 5 => 3 4 5 => 3 4 5
 * C 6 7 8      7 8     7 8
 * This follows the H-Gui convention, not the 'standard' convention
 */

/* neighbours are laid out in this pattern:
 *      6  12   7
 *   17   0   1  13
 * 11   5   X   2   8
 *   16   4   3  14
 *     10  15   9
 */
const MoveScore neighbours[18] = {
	MoveScore(-1,-1, 3), MoveScore(0,-1, 3), MoveScore(1, 0, 3), MoveScore(1, 1, 3), MoveScore( 0, 1, 3), MoveScore(-1, 0, 3), //direct neighbours, clockwise
	MoveScore(-2,-2, 1), MoveScore(0,-2, 1), MoveScore(2, 0, 1), MoveScore(2, 2, 1), MoveScore( 0, 2, 1), MoveScore(-2, 0, 1), //corners of ring 2, easy to block
	MoveScore(-1,-2, 2), MoveScore(1,-1, 2), MoveScore(2, 1, 2), MoveScore(1, 2, 2), MoveScore(-1, 1, 2), MoveScore(-2,-1, 2), //sides of ring 2, virtual connections
	};

static MoveValid * staticneighbourlist[11] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}; //one per boardsize


class Board{
public:

	static const int default_size = 8;
	static const int min_size = 3;
	static const int max_size = 10;
	static const int max_vecsize = 19*19;

	static const int pattern_cells = 18;
	typedef uint64_t Pattern;

	struct Cell {
		Side    piece;   //who controls this cell, 0 for none, 1,2 for players
		uint8_t size;    //size of this group of cells
mutable uint16_t parent; //parent for this group of cells
		uint8_t corner;  //which corners are this group connected to
		uint8_t edge;    //which edges are this group connected to
mutable uint8_t mark;    //when doing a ring search, has this position been seen?
		uint8_t perm;    //is this a permanent piece or a randomly placed piece?
		Pattern pattern; //the pattern of pieces for neighbours, but from their perspective. Rotate 180 for my perpective

		Cell() : piece(Side::NONE), size(0), parent(0), corner(0), edge(0), mark(0), perm(0), pattern(0) { }
		Cell(Side p, unsigned int a, unsigned int s, unsigned int c, unsigned int e, Pattern t) :
			piece(p), size(s), parent(a), corner(c), edge(e), mark(0), perm(0), pattern(t) { }

		int numcorners() const { return BitsSetTable256[corner]; }
		int numedges()   const { return BitsSetTable256[edge];   }

		std::string to_s(int i) const {
			return "Cell " + to_str(i) +": "
				"piece: " + to_str(piece.to_i())+
				", size: " + to_str((int)size) +
				", parent: " + to_str((int)parent) +
				", corner: " + to_str((int)corner) + "/" + to_str(numcorners()) +
				", edge: " + to_str((int)edge) + "/" + to_str(numedges()) +
				", perm: " + to_str((int)perm) +
				", pattern: " + to_str((int)pattern);
		}
	};

	class MoveIterator { //only returns valid moves...
		const Board & board;
		int lineend;
		MoveValid move;
		bool unique;
		HashSet hashes;
	public:
		MoveIterator(const Board & b, bool Unique) : board(b), lineend(0), move(Move(M_SWAP), -1), unique(Unique) {
			if(board.outcome >= Outcome::DRAW){
				move = MoveValid(0, board.get_size_d(), -1); //already done
			} else {
				if(unique)
					hashes.init(board.movesremain());
				++(*this); //find the first valid move
			}
		}

		const MoveValid & operator * ()  const { return move; }
		const MoveValid * operator -> () const { return & move; }
		bool done() const { return (move.y >= board.get_size_d()); }
		bool operator == (const Board::MoveIterator & rhs) const { return (move == rhs.move); }
		bool operator != (const Board::MoveIterator & rhs) const { return (move != rhs.move); }
		MoveIterator & operator ++ (){ //prefix form
			while(true){
				do{
					move.x++;
					move.xy++;

					if(move.x >= lineend){
						move.y++;
						if(move.y >= board.get_size_d()){ //done
							move.xy = -1;
							return *this;
						}
						move.x = board.linestart(move.y);
						move.xy = board.xy(move.x, move.y);
						lineend = board.lineend(move.y);
					}
				}while(!board.valid_move_fast(move));

				if(unique){
					uint64_t h = board.test_hash(move, board.toplay());
					if(!hashes.add(h))
						continue;
				}
				break;
			}
			return *this;
		}
	};

private:
	char size; //the length of one side of the hexagon
	char sizem1; //size - 1
	char size_d; //diameter of the board = size*2-1

	short num_cells;
	short nummoves;
	short unique_depth; //update and test rotations/symmetry with less than this many pieces on the board
	Move last;
	Side toPlay;
	Outcome outcome;
	char wintype; //0 no win, 1 = edge, 2 = corner, 3 = ring

	std::vector<Cell> cells;
	Zobrist<12> hash;
	const MoveValid * neighbourlist;

public:
	bool check_rings; // whether to look for rings at all
	int perm_rings;   // how many permanent stones are needed for a ring to count

	Board(){
		size = 0;
	}

	Board(int s){
		size = s;
		sizem1 = s - 1;
		size_d = s*2-1;
		last = M_NONE;
		nummoves = 0;
		unique_depth = 5;
		toPlay = Side::P1;
		outcome = Outcome::UNKNOWN;
		wintype = 0;
		check_rings = true;
		perm_rings = 0;
		neighbourlist = get_neighbour_list();
		num_cells = vecsize() - size*sizem1;

		cells.resize(vecsize());

		for(int y = 0; y < size_d; y++){
			for(int x = 0; x < size_d; x++){
				int posxy = xy(x, y);
				Pattern p = 0, j = 3;
				for(const MoveValid * i = nb_begin(posxy), *e = nb_end_big_hood(i); i < e; i++){
					if(!i->onboard())
						p |= j;
					j <<= 2;
				}
				Side s = (onboard(x, y) ? Side::NONE : Side::UNDEF);
				cells[posxy] = Cell(s, posxy, 1, (1 << iscorner(x, y)), (1 << isedge(x, y)), pattern_reverse(p));
			}
		}
	}

	int memsize() const { return sizeof(Board) + sizeof(Cell)*vecsize(); }

	int get_size_d() const { return size_d; }
	int get_size() const{ return size; }

	int vecsize() const { return size_d*size_d; }
	int numcells() const { return num_cells; }

	int num_moves() const { return nummoves; }
	int movesremain() const { return (won() >= Outcome::DRAW ? 0 : num_cells - nummoves); }

	int xy(int x, int y)   const { return   y*size_d +   x; }
	int xy(const Move & m) const { return m.y*size_d + m.x; }
	int xy(const MoveValid & m) const { return m.xy; }

	int xyc(int x, int y)   const { return xy(  x + sizem1,   y + sizem1); }
	int xyc(const Move & m) const { return xy(m.x + sizem1, m.y + sizem1); }

	MoveValid yx(int i) const { return MoveValid(i % size, i / size, i); }

	const Cell * cell(int i)          const { return & cells[i]; }
	const Cell * cell(int x, int y)   const { return cell(xy(x,y)); }
	const Cell * cell(const Move & m) const { return cell(xy(m)); }
	const Cell * cell(const MoveValid & m) const { return cell(m.xy); }


	//assumes valid x,y
	Side get(int i)          const { return cells[i].piece; }
	Side get(int x, int y)   const { return get(xy(x, y)); }
	Side get(const Move & m) const { return get(xy(m)); }
	Side get(const MoveValid & m) const { return get(m.xy); }

	Side geton(const MoveValid & m) const { return (m.onboard() ? get(m.xy) : Side::UNDEF); }

	int local(const Move & m, Side turn) const { return local(xy(m), turn); }
	int local(int i,          Side turn) const {
		Pattern p = pattern(i);
		Pattern x = ((p & 0xAAAAAAAAAull) >> 1) ^ (p & 0x555555555ull); // p1 is now when p1 or p2 but not both (ie off the board)
		p = x & (turn == Side::P1 ? p : p >> 1); // now just the selected player
		return (p & 0x000000FFF ? 3 : 0) |
		       (p & 0x000FFF000 ? 2 : 0) |
		       (p & 0xFFF000000 ? 1 : 0);
	}


	//assumes x, y are in array bounds
	bool onboard_fast(int x, int y)   const { return (  y -   x < size) && (  x -   y < size); }
	bool onboard_fast(const Move & m) const { return (m.y - m.x < size) && (m.x - m.y < size); }
	//checks array bounds too
	bool onboard(int x, int y)  const { return (  x >= 0 &&   y >= 0 &&   x < size_d &&   y < size_d && onboard_fast(x, y) ); }
	bool onboard(const Move & m)const { return (m.x >= 0 && m.y >= 0 && m.x < size_d && m.y < size_d && onboard_fast(m) ); }
	bool onboard(const MoveValid & m) const { return m.onboard(); }

	//assumes x, y are in bounds and the game isn't already finished
	bool valid_move_fast(int i)               const { return get(i) == Side::NONE; }
	bool valid_move_fast(int x, int y)        const { return valid_move_fast(xy(x, y)); }
	bool valid_move_fast(const Move & m)      const { return valid_move_fast(xy(m)); }
	bool valid_move_fast(const MoveValid & m) const { return valid_move_fast(m.xy); }
	//checks array bounds too
	bool valid_move(int x, int y)        const { return (outcome < Outcome::DRAW && onboard(x, y) && valid_move_fast(x, y)); }
	bool valid_move(const Move & m)      const { return (outcome < Outcome::DRAW && onboard(m)    && valid_move_fast(m)); }
	bool valid_move(const MoveValid & m) const { return (outcome < Outcome::DRAW && m.onboard()   && valid_move_fast(m)); }

	//iterator through neighbours of a position
	const MoveValid * nb_begin(int x, int y)   const { return nb_begin(xy(x, y)); }
	const MoveValid * nb_begin(const Move & m) const { return nb_begin(xy(m)); }
	const MoveValid * nb_begin(int i)          const { return &neighbourlist[i*18]; }

	const MoveValid * nb_end(int x, int y)   const { return nb_end(xy(x, y)); }
	const MoveValid * nb_end(const Move & m) const { return nb_end(xy(m)); }
	const MoveValid * nb_end(int i)          const { return nb_end(nb_begin(i)); }
	const MoveValid * nb_end(const MoveValid * m) const { return m + 6; }
	const MoveValid * nb_end_small_hood(const MoveValid * m) const { return m + 12; }
	const MoveValid * nb_end_big_hood(const MoveValid * m) const { return m + 18; }

	int iscorner(int x, int y) const {
		if(!onboard(x,y))
			return -1;

		int m = sizem1, e = size_d-1;

		if(x == 0 && y == 0) return 0;
		if(x == m && y == 0) return 1;
		if(x == e && y == m) return 2;
		if(x == e && y == e) return 3;
		if(x == m && y == e) return 4;
		if(x == 0 && y == m) return 5;

		return -1;
	}

	int isedge(int x, int y) const {
		if(!onboard(x,y))
			return -1;

		int m = sizem1, e = size_d-1;

		if(y   == 0 && x != 0 && x != m) return 0;
		if(x-y == m && x != m && x != e) return 1;
		if(x   == e && y != m && y != e) return 2;
		if(y   == e && x != e && x != m) return 3;
		if(y-x == m && x != m && x != 0) return 4;
		if(x   == 0 && y != m && y != 0) return 5;

		return -1;
	}

	MoveValid * get_neighbour_list(){
		if(!staticneighbourlist[(int)size]){
			MoveValid * list = new MoveValid[vecsize()*18];
			MoveValid * a = list;
			for(int y = 0; y < size_d; y++){
				for(int x = 0; x < size_d; x++){
					Move pos(x,y);

					for(int i = 0; i < 18; i++){
						Move loc = pos + neighbours[i];
						*a = MoveValid(loc, (onboard(loc) ? xy(loc) : -1) );
						++a;
					}
				}
			}

			staticneighbourlist[(int)size] = list;
		}

		return staticneighbourlist[(int)size];
	}


	int linestart(int y) const { return (y < size ? 0 : y - sizem1); }
	int lineend(int y)   const { return (y < size ? size + y : size_d); }
	int linelen(int y)   const { return size_d - abs(sizem1 - y); }

	std::string to_s(bool color) const {
		using std::string;
		string white = "O",
		       black = "@",
		       empty = ".",
		       coord = "",
		       reset = "";
		if(color){
			string esc = "\033";
			reset = esc + "[0m";
			coord = esc + "[1;37m";
			empty = reset + ".";
			white = esc + "[1;33m" + "@"; //yellow
			black = esc + "[1;34m" + "@"; //blue
		}

		string s;
		s += string(size + 3, ' ');
		for(int i = 0; i < size; i++)
			s += " " + coord + to_str(i+1);
		s += "\n";

		for(int y = 0; y < size_d; y++){
			s += string(abs(sizem1 - y) + 2, ' ');
			s += coord + char('A' + y);
			int end = lineend(y);
			for(int x = linestart(y); x < end; x++){
				s += (last == Move(x, y)   ? coord + "[" :
				      last == Move(x-1, y) ? coord + "]" : " ");
				Side p = get(x, y);
				if(     p == Side::NONE) s += empty;
				else if(p == Side::P1)   s += white;
				else if(p == Side::P2)   s += black;
				else                     s += "?";
			}
			s += (last == Move(end-1, y) ? coord + "]" : " ");
			if(y < sizem1)
				s += coord + to_str(size + y + 1);
			s += '\n';
		}

		s += reset;
		return s;
	}

	friend std::ostream& operator<< (std::ostream &out, const Board & b) { return out << b.to_s(true); }
	void print(bool color = true) const {
		printf("%s", to_s(color).c_str());
	}

	Outcome won() const {
		return outcome;
	}

	char getwintype() const { return wintype; }

	Side toplay() const {
		return toPlay;
	}

	MoveIterator moveit(bool unique = false) const {
		return MoveIterator(*this, (unique ? nummoves <= unique_depth : false));
	}

	void set(const Move & m, bool perm = true){
		last = m;
		Cell * cell = & cells[xy(m)];
		cell->piece = toPlay;
		cell->perm = perm;
		nummoves++;
		update_hash(m, toPlay); //depends on nummoves
		toPlay = ~toPlay;
	}

	void unset(const Move & m){ //break win checks, but is a poor mans undo if all you care about is the hash
		toPlay = ~toPlay;
		update_hash(m, toPlay);
		nummoves--;
		Cell * cell = & cells[xy(m)];
		cell->piece = Side::NONE;
		cell->perm = 0;
	}

	int find_group(const MoveValid & m) const { return find_group(m.xy); }
	int find_group(const Move & m) const { return find_group(xy(m)); }
	int find_group(int x, int y)   const { return find_group(xy(x, y)); }
	int find_group(unsigned int i) const {
		unsigned int p = cells[i].parent;
		if(p != i){
			do{
				p = cells[p].parent;
			}while(p != cells[p].parent);
			cells[i].parent = p; //do path compression, but only the current one, not all, to avoid recursion
		}
		return p;
	}

	//join the groups of two positions, propagating group size, and edge/corner connections
	//returns true if they're already the same group, false if they are now joined
	bool join_groups(const Move & a, const Move & b) { return join_groups(xy(a), xy(b)); }
	bool join_groups(int x1, int y1, int x2, int y2) { return join_groups(xy(x1, y1), xy(x2, y2)); }
	bool join_groups(int i, int j){
		i = find_group(i);
		j = find_group(j);

		if(i == j)
			return true;

		if(cells[i].size < cells[j].size) //force i's subtree to be bigger
			std::swap(i, j);

		cells[j].parent = i;
		cells[i].size   += cells[j].size;
		cells[i].corner |= cells[j].corner;
		cells[i].edge   |= cells[j].edge;

		return false;
	}

	Cell test_cell(const Move & pos) const {
		Side turn = toplay();
		int posxy = xy(pos);

		Cell testcell = cells[find_group(pos)];
		for(const MoveValid * i = nb_begin(posxy), *e = nb_end(i); i < e; i++){
			if(i->onboard() && turn == get(i->xy)){
				const Cell * g = & cells[find_group(i->xy)];
				testcell.corner |= g->corner;
				testcell.edge   |= g->edge;
				testcell.size   += g->size; //not quite accurate if it's joining the same group twice
				i++; //skip the next one
			}
		}
		return testcell;
	}

	int test_connectivity(const Move & pos) const {
		Cell testcell = test_cell(pos);
		return testcell.numcorners() + testcell.numedges();
	}

	int test_size(const Move & pos) const {
		Cell testcell = test_cell(pos);
		return testcell.size;
	}

	//check if a position is encirclable by a given player
	//false if it or one of its neighbours are the opponent's and connected to an edge or corner
	bool encirclable(const Move pos, Side player) const {
		Side otherplayer = ~player;
		int posxy = xy(pos);

		const Cell * g = & cells[find_group(posxy)];
		if(g->piece == otherplayer && (g->edge || g->corner))
			return false;

		for(const MoveValid * i = nb_begin(posxy), *e = nb_end(i); i < e; i++){
			if(!i->onboard())
				return false;

			const Cell * g = & cells[find_group(i->xy)];
			if(g->piece == otherplayer && (g->edge || g->corner))
				return false;
		}
		return true;
	}

	// do a depth first search for a ring
	bool checkring_df(const Move & pos, const Side turn) const {
		const Cell * start = & cells[xy(pos)];
		start->mark = 1;
		bool success = false;
		for(int i = 0; i < 4; i++){ //4 instead of 6 since any ring must have its first endpoint in the first 4
			Move loc = pos + neighbours[i];

			if(!onboard(loc))
				continue;

			const Cell * g = & cells[xy(loc)];

			if(turn != g->piece)
				continue;

			g->mark = 1;
			success = followring(loc, i, turn, (perm_rings - g->perm));
			g->mark = 0;

			if(success)
				break;
		}
		start->mark = 0;
		return success;
	}
	// only take the 3 directions that are valid in a ring
	// the backwards directions are either invalid or not part of the shortest loop
	bool followring(const Move & cur, const int & dir, const Side & turn, const int & permsneeded) const {
		for(int i = 5; i <= 7; i++){
			int nd = (dir + i) % 6;
			Move next = cur + neighbours[nd];

			if(!onboard(next))
				continue;

			const Cell * g = & cells[xy(next)];

			if(g->mark)
				return (permsneeded <= 0);

			if(turn != g->piece)
				continue;

			g->mark = 1;
			bool success = followring(next, nd, turn, (permsneeded - g->perm));
			g->mark = 0;

			if(success)
				return true;
		}
		return false;
	}

	// do an O(1) ring check
	// must be done before placing the stone and joining it with the neighbouring groups
	bool checkring_o1(const Move & pos, const Side turn) const {
		static const unsigned char ringdata[64][10] = {
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000000
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000001
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000010
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000011
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000100
			{1, 3, 5, 0, 0, 0, 0, 0, 0, 0}, //000101
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //000110
			{3,10,16,15, 0, 0, 0, 0, 0, 0}, //000111
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //001000
			{1, 2, 5, 0, 0, 0, 0, 0, 0, 0}, //001001
			{1, 2, 4, 0, 0, 0, 0, 0, 0, 0}, //001010
			{1, 2, 4, 0, 0, 0, 0, 0, 0, 0}, //001011
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //001100
			{1, 2, 5, 0, 0, 0, 0, 0, 0, 0}, //001101
			{3, 9,15,14, 0, 0, 0, 0, 0, 0}, //001110
			{4,10,16,15, 9,14,15, 0, 0, 0}, //001111
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //010000
			{1, 1, 5, 0, 0, 0, 0, 0, 0, 0}, //010001
			{1, 1, 4, 0, 0, 0, 0, 0, 0, 0}, //010010
			{1, 1, 4, 0, 0, 0, 0, 0, 0, 0}, //010011
			{1, 1, 3, 0, 0, 0, 0, 0, 0, 0}, //010100
			{2, 1, 3, 5, 0, 0, 0, 0, 0, 0}, //010101
			{1, 1, 3, 0, 0, 0, 0, 0, 0, 0}, //010110
			{7,10,16,15, 1, 3, 0, 0, 0, 0}, //010111
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //011000
			{1, 1, 5, 0, 0, 0, 0, 0, 0, 0}, //011001
			{1, 1, 4, 0, 0, 0, 0, 0, 0, 0}, //011010
			{1, 1, 4, 0, 0, 0, 0, 0, 0, 0}, //011011
			{3, 8,14,13, 0, 0, 0, 0, 0, 0}, //011100
			{7, 8,14,13, 1, 5, 0, 0, 0, 0}, //011101
			{4, 9,15,14, 8,13,14, 0, 0, 0}, //011110
			{5,10,16,15, 9,14,15, 8,14,13}, //011111
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //100000
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //100001
			{1, 0, 4, 0, 0, 0, 0, 0, 0, 0}, //100010
			{3,11,17,16, 0, 0, 0, 0, 0, 0}, //100011
			{1, 0, 3, 0, 0, 0, 0, 0, 0, 0}, //100100
			{1, 0, 3, 0, 0, 0, 0, 0, 0, 0}, //100101
			{1, 0, 3, 0, 0, 0, 0, 0, 0, 0}, //100110
			{4,11,17,16,10,15,16, 0, 0, 0}, //100111
			{1, 0, 2, 0, 0, 0, 0, 0, 0, 0}, //101000
			{1, 0, 2, 0, 0, 0, 0, 0, 0, 0}, //101001
			{2, 0, 2, 4, 0, 0, 0, 0, 0, 0}, //101010
			{7,11,17,16, 0, 2, 0, 0, 0, 0}, //101011
			{1, 0, 2, 0, 0, 0, 0, 0, 0, 0}, //101100
			{1, 0, 2, 0, 0, 0, 0, 0, 0, 0}, //101101
			{7, 9,15,14, 0, 2, 0, 0, 0, 0}, //101110
			{5,11,17,16,10,15,16, 9,15,14}, //101111
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //110000
			{3, 6,12,17, 0, 0, 0, 0, 0, 0}, //110001
			{1, 0, 4, 0, 0, 0, 0, 0, 0, 0}, //110010
			{4, 6,12,17,11,16,17, 0, 0, 0}, //110011
			{1, 0, 3, 0, 0, 0, 0, 0, 0, 0}, //110100
			{7, 6,12,17, 0, 3, 0, 0, 0, 0}, //110101
			{1, 0, 3, 0, 0, 0, 0, 0, 0, 0}, //110110
			{5, 6,12,17,11,16,17,10,16,15}, //110111
			{3, 7,13,12, 0, 0, 0, 0, 0, 0}, //111000
			{4, 7,13,12, 6,17,12, 0, 0, 0}, //111001
			{7, 7,13,12, 0, 4, 0, 0, 0, 0}, //111010
			{5, 7,13,12, 6,17,12,11,17,16}, //111011
			{4, 8,14,13, 7,12,13, 0, 0, 0}, //111100
			{5, 8,14,13, 7,12,13, 6,12,17}, //111101
			{5, 9,15,14, 8,13,14, 7,13,12}, //111110
			{6, 0, 0, 0, 0, 0, 0, 0, 0, 0}, //111111
		};

		int bitpattern = 0;
		const MoveValid * s = nb_begin(pos);
		for(const MoveValid * i = s, *e = nb_end(i); i < e; i++){
			bitpattern <<= 1;
			if(i->onboard() && turn == get(i->xy))
				bitpattern |= 1;
		}

		const unsigned char * d = ringdata[bitpattern];

		switch(d[0]){
			case 0: //no ring (000000, 000001, 000011)
				return false;

			case 1: //simple case (000101, 001101, 001011, 011011)
				return (find_group(s[d[1]]) == find_group(s[d[2]]));

			case 2:{ //3 non-neighbours (010101)
				int a = find_group(s[d[1]]), b = find_group(s[d[2]]), c = find_group(s[d[3]]);
				return (a == b || a == c || b == c);
			}

			case 7: //case 1 and 3 (010111)
				if(find_group(s[d[4]]) == find_group(s[d[5]]))
					return true;
				//fall through

			case 3: // 3 neighbours (000111)
				return checkring_back(s[d[1]], s[d[2]], s[d[3]], turn);

			case 4: // 4 neighbours (001111)
				return checkring_back(s[d[1]], s[d[2]], s[d[3]], turn) ||
				       checkring_back(s[d[4]], s[d[5]], s[d[6]], turn);

			case 5: // 5 neighbours (011111)
				return checkring_back(s[d[1]], s[d[2]], s[d[3]], turn) ||
				       checkring_back(s[d[4]], s[d[5]], s[d[6]], turn) ||
				       checkring_back(s[d[7]], s[d[8]], s[d[9]], turn);

			case 6: // 6 neighbours (111111)
				return true; //a ring around this position? how'd that happen

			default:
				return false;
		}
	}
	//checks for 3 more stones, a should be the corner
	bool checkring_back(const MoveValid & a, const MoveValid & b, const MoveValid & c, Side turn) const {
		return (a.onboard() && get(a) == turn && get(b) == turn && get(c) == turn);
	}

	hash_t gethash() const {
		return (nummoves > unique_depth ? hash.get(0) : hash.get());
	}

	std::string hashstr() const {
		static const char hexlookup[] = "0123456789abcdef";
		char buf[19] = "0x";
		hash_t val = gethash();
		for(int i = 15; i >= 0; i--){
			buf[i+2] = hexlookup[val & 15];
			val >>= 4;
		}
		buf[18] = '\0';
		return (char *)buf;
	}

	void update_hash(const Move & pos, Side side){
		int turn = side.to_i();
		if(nummoves > unique_depth){ //simple update, no rotations/symmetry
			hash.update(0, 3*xy(pos) + turn);
			return;
		}

		//mirror is simply flip x,y
		int x = pos.x - sizem1,
		    y = pos.y - sizem1,
		    z = y - x;

//x,y; y,z; z,-x; -x,-y; -y,-z; -z,x
//y,x; z,y; -x,z; -y,-x; -z,-y; x,-z

		hash.update(0,  3*xyc( x,  y) + turn);
		hash.update(1,  3*xyc( y,  z) + turn);
		hash.update(2,  3*xyc( z, -x) + turn);
		hash.update(3,  3*xyc(-x, -y) + turn);
		hash.update(4,  3*xyc(-y, -z) + turn);
		hash.update(5,  3*xyc(-z,  x) + turn);
		hash.update(6,  3*xyc( y,  x) + turn);
		hash.update(7,  3*xyc( z,  y) + turn);
		hash.update(8,  3*xyc(-x,  z) + turn);
		hash.update(9,  3*xyc(-y, -x) + turn);
		hash.update(10, 3*xyc(-z, -y) + turn);
		hash.update(11, 3*xyc( x, -z) + turn);
	}

	hash_t test_hash(const Move & pos) const {
		return test_hash(pos, toplay());
	}

	hash_t test_hash(const Move & pos, Side side) const {
		int turn = side.to_i();
		if(nummoves >= unique_depth) //simple test, no rotations/symmetry
			return hash.test(0, 3*xy(pos) + turn);

		int x = pos.x - sizem1,
		    y = pos.y - sizem1,
		    z = y - x;

		hash_t m = hash.test(0,  3*xyc( x,  y) + turn);
		m = std::min(m, hash.test(1,  3*xyc( y,  z) + turn));
		m = std::min(m, hash.test(2,  3*xyc( z, -x) + turn));
		m = std::min(m, hash.test(3,  3*xyc(-x, -y) + turn));
		m = std::min(m, hash.test(4,  3*xyc(-y, -z) + turn));
		m = std::min(m, hash.test(5,  3*xyc(-z,  x) + turn));
		m = std::min(m, hash.test(6,  3*xyc( y,  x) + turn));
		m = std::min(m, hash.test(7,  3*xyc( z,  y) + turn));
		m = std::min(m, hash.test(8,  3*xyc(-x,  z) + turn));
		m = std::min(m, hash.test(9,  3*xyc(-y, -x) + turn));
		m = std::min(m, hash.test(10, 3*xyc(-z, -y) + turn));
		m = std::min(m, hash.test(11, 3*xyc( x, -z) + turn));
		return m;
	}

	Pattern sympattern(const MoveValid & pos) const { return sympattern(pos.xy); }
	Pattern sympattern(const Move & pos)      const { return sympattern(xy(pos)); }
	Pattern sympattern(int posxy)             const { return pattern_symmetry(pattern(posxy)); }

	Pattern pattern(const MoveValid & pos) const { return pattern(pos.xy); }
	Pattern pattern(const Move & pos)      const { return pattern(xy(pos)); }
	Pattern pattern(int posxy)             const {
		// this is from the opposite perspective
		// so rotate into this move's perspective
		return pattern_reverse(cells[posxy].pattern);
	}

	Pattern pattern_medium(const MoveValid & pos) const { return pattern_medium(pos.xy); }
	Pattern pattern_medium(const Move & pos)      const { return pattern_medium(xy(pos)); }
	Pattern pattern_medium(int posxy)             const {
		return pattern(posxy) & ((1ull << 24) - 1);
	}

	Pattern pattern_small(const MoveValid & pos) const { return pattern_small(pos.xy); }
	Pattern pattern_small(const Move & pos)      const { return pattern_small(xy(pos)); }
	Pattern pattern_small(int posxy)             const {
		return pattern(posxy) & ((1ull << 12) - 1);
	}

	static Pattern pattern_reverse(Pattern p) { // switch perspectives (position out, or position in)
		return (((p & 0x03F03F03Full) << 6) | ((p & 0xFC0FC0FC0ull) >> 6));
	}

	static Pattern pattern_invert(Pattern p){ //switch players
		return ((p & 0xAAAAAAAAAull) >> 1) | ((p & 0x555555555ull) << 1);
	}
	static Pattern pattern_rotate(Pattern p){
		return (((p & 0x003003003ull) << 10) | ((p & 0xFFCFFCFFCull) >> 2));
	}
	static Pattern pattern_mirror(Pattern p){
		// HGFEDC BA9876 543210 -> DEFGHC 6789AB 123450
		return ((p & (3ull <<  6))      ) | ((p & (3ull <<  0))     ) | // 0,3 stay in place
		       ((p & (3ull << 10)) >>  8) | ((p & (3ull <<  2)) << 8) | // 1,5 swap
		       ((p & (3ull <<  8)) >>  4) | ((p & (3ull <<  4)) << 4) | // 2,4 swap
		       ((p & (3ull << 22)) >> 10) | ((p & (3ull << 12)) <<10) | // 6,B swap
		       ((p & (3ull << 20)) >>  6) | ((p & (3ull << 14)) << 6) | // 7,A swap
		       ((p & (3ull << 18)) >>  2) | ((p & (3ull << 16)) << 2) | // 8,9 swap
		       ((p & (3ull << 30))      ) | ((p & (3ull << 24))     ) | // F,C stay in place
		       ((p & (3ull << 34)) >>  8) | ((p & (3ull << 26)) << 8) | // H,D swap
		       ((p & (3ull << 32)) >>  4) | ((p & (3ull << 28)) << 4);  // G,E swap
	}
	static Pattern pattern_symmetry(Pattern p){ //takes a pattern and returns the representative version
		Pattern m = p;                 //012345
		m = std::min(m, (p = pattern_rotate(p)));//501234
		m = std::min(m, (p = pattern_rotate(p)));//450123
		m = std::min(m, (p = pattern_rotate(p)));//345012
		m = std::min(m, (p = pattern_rotate(p)));//234501
		m = std::min(m, (p = pattern_rotate(p)));//123450
		m = std::min(m, (p = pattern_mirror(pattern_rotate(p))));//012345 -> 054321
		m = std::min(m, (p = pattern_rotate(p)));//105432
		m = std::min(m, (p = pattern_rotate(p)));//210543
		m = std::min(m, (p = pattern_rotate(p)));//321054
		m = std::min(m, (p = pattern_rotate(p)));//432105
		m = std::min(m, (p = pattern_rotate(p)));//543210
		return m;
	}

	bool move(const Move & pos, bool checkwin = true, bool permanent = true){
		return move(MoveValid(pos, xy(pos)), checkwin, permanent);
	}
	bool move(const MoveValid & pos, bool checkwin = true, bool permanent = true){
		assert(outcome < Outcome::DRAW);

		if(!valid_move(pos))
			return false;

		Side turn = toplay();
		set(pos, permanent);

		// update the nearby patterns
		Pattern p = turn.to_i();
		for(const MoveValid * i = nb_begin(pos.xy), *e = nb_end_big_hood(i); i < e; i++){
			if(i->onboard()){
				cells[i->xy].pattern |= p;
			}
			p <<= 2;
		}

		// join the groups for win detection
		bool alreadyjoined = false; //useful for finding rings
		for(const MoveValid * i = nb_begin(pos.xy), *e = nb_end(i); i < e; i++){
			if(i->onboard() && turn == get(i->xy)){
				alreadyjoined |= join_groups(pos.xy, i->xy);
				i++; //skip the next one. If it is the same group,
					 //it is already connected and forms a corner, which we can ignore
			}
		}

		if(checkwin){
			Cell * g = & cells[find_group(pos.xy)];
			if(g->numedges() >= 3){
				outcome = +turn;
				wintype = 1;
			}else if(g->numcorners() >= 2){
				outcome = +turn;
				wintype = 2;
			}else if(check_rings && alreadyjoined && g->size >= 6 && checkring_df(pos, turn)){
				outcome = +turn;
				wintype = 3;
			}else if(nummoves == num_cells){
				outcome = Outcome::DRAW;
			}
		}
		return true;
	}

	bool test_local(const Move & pos, Side turn) const {
		return (local(pos, turn) == 3);
	}

	//test if making this move would win, but don't actually make the move
	Outcome test_outcome(const Move & pos) const { return test_outcome(pos, toplay()); }
	Outcome test_outcome(const Move & pos, Side turn) const {
		if(test_local(pos, turn)){
			int posxy = xy(pos);
			Cell testcell = cells[find_group(posxy)];
			int numgroups = 0;
			for(const MoveValid * i = nb_begin(posxy), *e = nb_end(i); i < e; i++){
				if(i->onboard() && turn == get(i->xy)){
					const Cell * g = & cells[find_group(i->xy)];
					testcell.corner |= g->corner;
					testcell.edge   |= g->edge;
					testcell.size   += g->size;
					i++; //skip the next one
					numgroups++;
				}
			}

			if(testcell.numcorners() >= 2 || testcell.numedges() >= 3 || (check_rings && numgroups >= 2 && testcell.size >= 6 && checkring_o1(pos, turn)))
				return +turn;
		}

		if(nummoves+1 == num_cells)
			return Outcome::DRAW;

		return Outcome::UNKNOWN;
	}
};

}; // namespace Havannah
}; // namespace Morat
