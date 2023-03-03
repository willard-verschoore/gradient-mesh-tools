#include "gmt/slot-map.hpp"

#include <catch2/catch.hpp>
#include <string>

using namespace gmt;

TEST_CASE("Insert/removal operations affect size", "[slot-map]")
{
  SlotMap<std::string> strings;
  REQUIRE(strings.size() == 0);
  auto first = strings.add("first item");
  REQUIRE(strings.size() == 1);
  auto second = strings.add("second item");
  REQUIRE(strings.size() == 2);
  strings.remove(second);
  REQUIRE(strings.size() == 1);
  strings.remove(first);
  REQUIRE(strings.size() == 0);
}

TEST_CASE("Accessing items through handles", "[slot-map]")
{
  SlotMap<std::string> strings;
  auto first = strings.add("first item");
  auto second = strings.add("second item");

  REQUIRE(*strings.find(first) == "first item");
  REQUIRE(*strings.find(second) == "second item");
}

TEST_CASE("Removing an item does not affect others", "[slot-map]")
{
  SlotMap<std::string> strings;
  auto first = strings.add("first item");
  auto second = strings.add("second item");

  REQUIRE(*strings.find(first) == "first item");
  REQUIRE(*strings.find(second) == "second item");

  strings.remove(first);
  REQUIRE(strings.find(first) == strings.end());
  REQUIRE(*strings.find(second) == "second item");
}
