#include "grid.hpp"
#include <gtest/gtest.h>

using Grid = FsGrid<float, 1>;

TEST(FsGridTest, xyzToLinear) {
   ASSERT_EQ(0, Grid::xyzToLinear(-1, -1, -1));
   ASSERT_EQ(1, Grid::xyzToLinear(-1, -1, 0));
   ASSERT_EQ(2, Grid::xyzToLinear(-1, -1, 1));
   ASSERT_EQ(3, Grid::xyzToLinear(-1, 0, -1));
   ASSERT_EQ(4, Grid::xyzToLinear(-1, 0, 0));
   ASSERT_EQ(5, Grid::xyzToLinear(-1, 0, 1));
   ASSERT_EQ(6, Grid::xyzToLinear(-1, 1, -1));
   ASSERT_EQ(7, Grid::xyzToLinear(-1, 1, 0));
   ASSERT_EQ(8, Grid::xyzToLinear(-1, 1, 1));
   ASSERT_EQ(9, Grid::xyzToLinear(0, -1, -1));
   ASSERT_EQ(10, Grid::xyzToLinear(0, -1, 0));
   ASSERT_EQ(11, Grid::xyzToLinear(0, -1, 1));
   ASSERT_EQ(12, Grid::xyzToLinear(0, 0, -1));
   ASSERT_EQ(13, Grid::xyzToLinear(0, 0, 0));
   ASSERT_EQ(14, Grid::xyzToLinear(0, 0, 1));
   ASSERT_EQ(15, Grid::xyzToLinear(0, 1, -1));
   ASSERT_EQ(16, Grid::xyzToLinear(0, 1, 0));
   ASSERT_EQ(17, Grid::xyzToLinear(0, 1, 1));
   ASSERT_EQ(18, Grid::xyzToLinear(1, -1, -1));
   ASSERT_EQ(19, Grid::xyzToLinear(1, -1, 0));
   ASSERT_EQ(20, Grid::xyzToLinear(1, -1, 1));
   ASSERT_EQ(21, Grid::xyzToLinear(1, 0, -1));
   ASSERT_EQ(22, Grid::xyzToLinear(1, 0, 0));
   ASSERT_EQ(23, Grid::xyzToLinear(1, 0, 1));
   ASSERT_EQ(24, Grid::xyzToLinear(1, 1, -1));
   ASSERT_EQ(25, Grid::xyzToLinear(1, 1, 0));
   ASSERT_EQ(26, Grid::xyzToLinear(1, 1, 1));
}

TEST(FsGridTest, linearToX) {
   ASSERT_EQ(-1, Grid::linearToX(0));
   ASSERT_EQ(-1, Grid::linearToX(1));
   ASSERT_EQ(-1, Grid::linearToX(2));
   ASSERT_EQ(-1, Grid::linearToX(3));
   ASSERT_EQ(-1, Grid::linearToX(4));
   ASSERT_EQ(-1, Grid::linearToX(5));
   ASSERT_EQ(-1, Grid::linearToX(6));
   ASSERT_EQ(-1, Grid::linearToX(7));
   ASSERT_EQ(-1, Grid::linearToX(8));
   ASSERT_EQ(0, Grid::linearToX(9));
   ASSERT_EQ(0, Grid::linearToX(10));
   ASSERT_EQ(0, Grid::linearToX(11));
   ASSERT_EQ(0, Grid::linearToX(12));
   ASSERT_EQ(0, Grid::linearToX(13));
   ASSERT_EQ(0, Grid::linearToX(14));
   ASSERT_EQ(0, Grid::linearToX(15));
   ASSERT_EQ(0, Grid::linearToX(16));
   ASSERT_EQ(0, Grid::linearToX(17));
   ASSERT_EQ(1, Grid::linearToX(18));
   ASSERT_EQ(1, Grid::linearToX(19));
   ASSERT_EQ(1, Grid::linearToX(20));
   ASSERT_EQ(1, Grid::linearToX(21));
   ASSERT_EQ(1, Grid::linearToX(22));
   ASSERT_EQ(1, Grid::linearToX(23));
   ASSERT_EQ(1, Grid::linearToX(24));
   ASSERT_EQ(1, Grid::linearToX(25));
   ASSERT_EQ(1, Grid::linearToX(26));
}

TEST(FsGridTest, linearToY) {
   ASSERT_EQ(-1, Grid::linearToY(0));
   ASSERT_EQ(-1, Grid::linearToY(1));
   ASSERT_EQ(-1, Grid::linearToY(2));
   ASSERT_EQ(0, Grid::linearToY(3));
   ASSERT_EQ(0, Grid::linearToY(4));
   ASSERT_EQ(0, Grid::linearToY(5));
   ASSERT_EQ(1, Grid::linearToY(6));
   ASSERT_EQ(1, Grid::linearToY(7));
   ASSERT_EQ(1, Grid::linearToY(8));
   ASSERT_EQ(-1, Grid::linearToY(9));
   ASSERT_EQ(-1, Grid::linearToY(10));
   ASSERT_EQ(-1, Grid::linearToY(11));
   ASSERT_EQ(0, Grid::linearToY(12));
   ASSERT_EQ(0, Grid::linearToY(13));
   ASSERT_EQ(0, Grid::linearToY(14));
   ASSERT_EQ(1, Grid::linearToY(15));
   ASSERT_EQ(1, Grid::linearToY(16));
   ASSERT_EQ(1, Grid::linearToY(17));
   ASSERT_EQ(-1, Grid::linearToY(18));
   ASSERT_EQ(-1, Grid::linearToY(19));
   ASSERT_EQ(-1, Grid::linearToY(20));
   ASSERT_EQ(0, Grid::linearToY(21));
   ASSERT_EQ(0, Grid::linearToY(22));
   ASSERT_EQ(0, Grid::linearToY(23));
   ASSERT_EQ(1, Grid::linearToY(24));
   ASSERT_EQ(1, Grid::linearToY(25));
   ASSERT_EQ(1, Grid::linearToY(26));
}

TEST(FsGridTest, linearToZ) {
   ASSERT_EQ(-1, Grid::linearToZ(0));
   ASSERT_EQ(0, Grid::linearToZ(1));
   ASSERT_EQ(1, Grid::linearToZ(2));
   ASSERT_EQ(-1, Grid::linearToZ(3));
   ASSERT_EQ(0, Grid::linearToZ(4));
   ASSERT_EQ(1, Grid::linearToZ(5));
   ASSERT_EQ(-1, Grid::linearToZ(6));
   ASSERT_EQ(0, Grid::linearToZ(7));
   ASSERT_EQ(1, Grid::linearToZ(8));
   ASSERT_EQ(-1, Grid::linearToZ(9));
   ASSERT_EQ(0, Grid::linearToZ(10));
   ASSERT_EQ(1, Grid::linearToZ(11));
   ASSERT_EQ(-1, Grid::linearToZ(12));
   ASSERT_EQ(0, Grid::linearToZ(13));
   ASSERT_EQ(1, Grid::linearToZ(14));
   ASSERT_EQ(-1, Grid::linearToZ(15));
   ASSERT_EQ(0, Grid::linearToZ(16));
   ASSERT_EQ(1, Grid::linearToZ(17));
   ASSERT_EQ(-1, Grid::linearToZ(18));
   ASSERT_EQ(0, Grid::linearToZ(19));
   ASSERT_EQ(1, Grid::linearToZ(20));
   ASSERT_EQ(-1, Grid::linearToZ(21));
   ASSERT_EQ(0, Grid::linearToZ(22));
   ASSERT_EQ(1, Grid::linearToZ(23));
   ASSERT_EQ(-1, Grid::linearToZ(24));
   ASSERT_EQ(0, Grid::linearToZ(25));
   ASSERT_EQ(1, Grid::linearToZ(26));
}

TEST(FsGridTest, xyzToLinearToxyz) {
   for (int32_t i = 0; i < 27; i++) {
      const auto x = Grid::linearToX(i);
      const auto y = Grid::linearToY(i);
      const auto z = Grid::linearToZ(i);
      ASSERT_EQ(i, Grid::xyzToLinear(x, y, z));
   }
}
