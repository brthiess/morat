
#include "../lib/catch.hpp"

#include "agentsolver.h"

using namespace Morat;
using namespace BridgIt;

TEST_CASE("BridgIt::AgentSolver::Node::to_s/from_s", "[hex][agentsolver]") {
	AgentSolver::Node n(Move("a1"));
	auto s = n.to_s();
	AgentSolver::Node k;
	REQUIRE(k.from_s(s));
	REQUIRE(n.to_s() == k.to_s());
}
