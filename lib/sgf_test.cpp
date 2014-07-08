
#include <sstream>

#include "catch.hpp"
#include "sgf.h"

#include "../havannah/move.h"

namespace Morat {

TEST_CASE("sgf simple", "[sgf]") {

	std::stringstream s;

	{ // write an sgf file
		SGFPrinter<Havannah::Move> sgf(s);
		sgf.game("havannah");
		sgf.size(5);
		sgf.end_root();
		sgf.move(Side::P1, Havannah::Move("a1"));
		sgf.comment("good");
		sgf.end();

		CHECK(s.str() == "(;FF[4]GM[havannah]SZ[5]\n"
		                 " ;W[a1]C[good]\n"
		                 ")\n");
	}

	{ // read one and get back what was written above
		SGFParser<Havannah::Move> sgf(s);
		REQUIRE(sgf.game() == "havannah");
		REQUIRE(sgf.size() == 5);
		REQUIRE(sgf.move() == Havannah::Move());
		REQUIRE(sgf.next_node());
		REQUIRE_FALSE(sgf.has_children());
		REQUIRE(sgf.move() == Havannah::Move("a1"));
		REQUIRE(sgf.comment() == "good");
		REQUIRE_FALSE(sgf.next_node());
		REQUIRE_FALSE(sgf.has_children());
	}
}

TEST_CASE("sgf write/read", "[sgf]") {
	std::stringstream s;

	{
		SGFPrinter<Havannah::Move> sgf(s);
		sgf.game("havannah");
		sgf.size(5);

		sgf.end_root();

		sgf.move(Side::P1, Havannah::Move("a1"));
		sgf.move(Side::P2, Havannah::Move("b2"));

		sgf.child_start();
			sgf.move(Side::P1, Havannah::Move("c1"));
			sgf.comment("c1");

			sgf.child_start();
				sgf.move(Side::P2, Havannah::Move("d1"));
				sgf.comment("d1");
			sgf.child_end();

		sgf.child_end();

		sgf.child_start();
			sgf.move(Side::P1, Havannah::Move("c2"));
			sgf.comment("c2");
		sgf.child_end();

		sgf.end();

		CHECK(s.str() == "(;FF[4]GM[havannah]SZ[5]\n"
		                 " ;W[a1];B[b2]\n"
		                 " (;W[c1]C[c1]\n"
		                 "  (;B[d1]C[d1])\n"
		                 " )\n"
		                 " (;W[c2]C[c2])\n"
		                 ")\n");
	}

	{ // read one and get back what was written above
		SGFParser<Havannah::Move> sgf(s);
		REQUIRE(sgf.game() == "havannah");
		REQUIRE(sgf.size() == 5);
		REQUIRE(sgf.move() == Havannah::Move());
		REQUIRE(sgf.next_node());
		REQUIRE(sgf.move() == Havannah::Move("a1"));
		REQUIRE(sgf.comment() == "");
		REQUIRE(sgf.next_node());
		REQUIRE(sgf.move() == Havannah::Move("b2"));
		REQUIRE_FALSE(sgf.next_node());
		REQUIRE(sgf.has_children());
		REQUIRE(sgf.next_child());
		REQUIRE(sgf.move() == Havannah::Move("c1"));
		REQUIRE(sgf.comment() == "c1");
		REQUIRE(sgf.has_children());
		REQUIRE(sgf.next_child());
		REQUIRE(sgf.move() == Havannah::Move("d1"));
		REQUIRE(sgf.comment() == "d1");
		REQUIRE_FALSE(sgf.has_children());
		REQUIRE_FALSE(sgf.next_child());
		REQUIRE(sgf.done_child());
		REQUIRE_FALSE(sgf.has_children());
		REQUIRE_FALSE(sgf.next_child());
		REQUIRE(sgf.done_child());
		REQUIRE(sgf.has_children());
		REQUIRE(sgf.next_child());
		REQUIRE(sgf.move() == Havannah::Move("c2"));
		REQUIRE(sgf.comment() == "c2");
		REQUIRE(sgf.done_child());
		REQUIRE(sgf.done_child());
	}
}

}; // namespace Morat
