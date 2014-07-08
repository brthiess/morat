
#include "../lib/catch.hpp"

#include "board.h"

using namespace Morat;
using namespace Chex;

void test_game(Board b, std::vector<std::string> moves, Outcome outcome) {
	REQUIRE(b.num_moves() == 0);
	Side side = Side::P1;
	for(auto s : moves) {
		Outcome expected = (s == moves.back() ? outcome : Outcome::UNKNOWN);
		Move move(s);
		CAPTURE(move);
		//CAPTURE(b);
		REQUIRE(b.valid_move(move));
		REQUIRE(b.toplay() == side);
		REQUIRE(b.test_outcome(move) == expected);
		REQUIRE(b.move(move));
		REQUIRE(b.won() == expected);
		side = ~side;
	}
}

TEST_CASE("Chex::Board", "[hex][board]") {
	Board b(3,6);

	SECTION("Basics") {
		REQUIRE(b.get_x_size() == 3);
		REQUIRE(b.get_y_size() == 6);		
		REQUIRE(b.movesremain() == 18);
	}

	SECTION("valid moves") {
		std::string valid[] = {
		"a1", "a2", "a3", 
		    "b1", "b2", "b3",
		       "c1", "c2", "c3",
		           "d1", "d2", "d3",
		               "e1", "e2", "e3",
		                   "f1", "f2", "f3", 
		};
		for(auto m : valid){
			REQUIRE(b.onboard(m));
			REQUIRE(b.valid_move(m));
		}
	}

	SECTION("invalid moves") {
		std::string invalid[] = {"a0", "a8", "a10", "b8", "c8", "e0", "e8", "f8", "f0", "h1", "f0"};
		for(auto m : invalid){
			REQUIRE_FALSE(b.onboard(m));
			REQUIRE_FALSE(b.valid_move(m));
		}
	}

	SECTION("duplicate moves") {
		Move m("a1");
		REQUIRE(b.valid_move(m));
		REQUIRE(b.move(m));
		REQUIRE_FALSE(b.valid_move(m));
		REQUIRE_FALSE(b.move(m));
	}

	SECTION("White Wins") {
		test_game(b, {      "a1", "b1", "a2", "b2", "a3", "b3",
							"c1", "d1", "c2", "d2", "c3", "d3",
							"e1", "f1", "e2", "f2", "e3", "f3"}, Outcome::P1);
	}

	SECTION("Black Wins") {
		test_game(b, {      "a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2", "e1", "e2", "f1", "f2"}, Outcome::P2);
		test_game(b, {"f1",	"a1", "a2", "b3", "b2", "c3", "c2", "d3", "d2", "e3", "e2", "f3"}, Outcome::P2);
	}

	SECTION("Unknown") {
		test_game(b, {      "b2", "f3" }, Outcome::UNKNOWN);
		test_game(b, {"d1", "b2", "f3", "b3"}, Outcome::UNKNOWN);
	}
}
