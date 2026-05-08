#include <gtest/gtest.h>
#include <hc/core/variable.hpp>

TEST(Variable, Constructor) {
  hc::Variable x("x");
  EXPECT_EQ(x.name, "x");
}

TEST(Variable, UnicodeNames) {
  hc::Variable alpha("α"), beta("β");
  EXPECT_EQ(alpha.name, "α");
  EXPECT_EQ(beta.name, "β");
}

TEST(MakeVariables, CreatesAndSorts) {
  auto vars = hc::make_variables({"x", "y", "z"});
  EXPECT_EQ(vars.size(), 3);
  EXPECT_EQ(vars[0].name, "x");
  EXPECT_EQ(vars[1].name, "y");
  EXPECT_EQ(vars[2].name, "z");
}

TEST(MergeSortedVars, DisjointSets) {
  auto v1 = hc::make_variables({"x", "z"});
  auto v2 = hc::make_variables({"y"});
  auto r = hc::merge_sorted_vars(v1, v2);
  EXPECT_EQ(r.vars, hc::make_variables({"x", "y", "z"}));
  EXPECT_EQ(r.map1, (std::vector<int>{0, 2}));
  EXPECT_EQ(r.map2, (std::vector<int>{1}));
}

TEST(MergeSortedVars, Overlap) {
  auto v1 = hc::make_variables({"x", "y"});
  auto v2 = hc::make_variables({"y", "z"});
  auto r = hc::merge_sorted_vars(v1, v2);
  EXPECT_EQ(r.vars, hc::make_variables({"x", "y", "z"}));
  EXPECT_EQ(r.map1, (std::vector<int>{0, 1}));
  EXPECT_EQ(r.map2, (std::vector<int>{1, 2}));
}

TEST(MergeSortedVars, IdenticalSets) {
  auto v1 = hc::make_variables({"x", "y"});
  auto v2 = hc::make_variables({"x", "y"});
  auto r = hc::merge_sorted_vars(v1, v2);
  EXPECT_EQ(r.vars, hc::make_variables({"x", "y"}));
  EXPECT_EQ(r.map1, (std::vector<int>{0, 1}));
  EXPECT_EQ(r.map2, (std::vector<int>{0, 1}));
}

TEST(MergeSortedVars, OneEmpty) {
  auto v1 = hc::make_variables({"x", "y"});
  auto v2 = std::vector<hc::Variable>{};
  auto r = hc::merge_sorted_vars(v1, v2);
  EXPECT_EQ(r.vars, hc::make_variables({"x", "y"}));
  EXPECT_EQ(r.map1, (std::vector<int>{0, 1}));
  EXPECT_TRUE(r.map2.empty());
}
