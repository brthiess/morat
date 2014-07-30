
#include "../lib/catch.hpp"

#include "agentmcts.h"

using namespace Morat;
using namespace BridgIt;

TEST_CASE("BridgIt::AgentMCTS::Node::to_s/from_s", "[hex][agentmcts]") {
	AgentMCTS::Node n(Move("a1"));
	auto s = n.to_s();
	AgentMCTS::Node k;
	REQUIRE(k.from_s(s));
	REQUIRE(n.to_s() == k.to_s());
}
