#include "grid.hpp"
#include <gtest/gtest.h>

TEST(FsGridTest, xyzToLinear) {
   ASSERT_EQ(0, fsgrid_detail::xyzToLinear(-1, -1, -1));
   ASSERT_EQ(1, fsgrid_detail::xyzToLinear(-1, -1, 0));
   ASSERT_EQ(2, fsgrid_detail::xyzToLinear(-1, -1, 1));
   ASSERT_EQ(3, fsgrid_detail::xyzToLinear(-1, 0, -1));
   ASSERT_EQ(4, fsgrid_detail::xyzToLinear(-1, 0, 0));
   ASSERT_EQ(5, fsgrid_detail::xyzToLinear(-1, 0, 1));
   ASSERT_EQ(6, fsgrid_detail::xyzToLinear(-1, 1, -1));
   ASSERT_EQ(7, fsgrid_detail::xyzToLinear(-1, 1, 0));
   ASSERT_EQ(8, fsgrid_detail::xyzToLinear(-1, 1, 1));
   ASSERT_EQ(9, fsgrid_detail::xyzToLinear(0, -1, -1));
   ASSERT_EQ(10, fsgrid_detail::xyzToLinear(0, -1, 0));
   ASSERT_EQ(11, fsgrid_detail::xyzToLinear(0, -1, 1));
   ASSERT_EQ(12, fsgrid_detail::xyzToLinear(0, 0, -1));
   ASSERT_EQ(13, fsgrid_detail::xyzToLinear(0, 0, 0));
   ASSERT_EQ(14, fsgrid_detail::xyzToLinear(0, 0, 1));
   ASSERT_EQ(15, fsgrid_detail::xyzToLinear(0, 1, -1));
   ASSERT_EQ(16, fsgrid_detail::xyzToLinear(0, 1, 0));
   ASSERT_EQ(17, fsgrid_detail::xyzToLinear(0, 1, 1));
   ASSERT_EQ(18, fsgrid_detail::xyzToLinear(1, -1, -1));
   ASSERT_EQ(19, fsgrid_detail::xyzToLinear(1, -1, 0));
   ASSERT_EQ(20, fsgrid_detail::xyzToLinear(1, -1, 1));
   ASSERT_EQ(21, fsgrid_detail::xyzToLinear(1, 0, -1));
   ASSERT_EQ(22, fsgrid_detail::xyzToLinear(1, 0, 0));
   ASSERT_EQ(23, fsgrid_detail::xyzToLinear(1, 0, 1));
   ASSERT_EQ(24, fsgrid_detail::xyzToLinear(1, 1, -1));
   ASSERT_EQ(25, fsgrid_detail::xyzToLinear(1, 1, 0));
   ASSERT_EQ(26, fsgrid_detail::xyzToLinear(1, 1, 1));
}

TEST(FsGridTest, linearToX) {
   ASSERT_EQ(-1, fsgrid_detail::linearToX(0));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(1));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(2));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(3));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(4));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(5));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(6));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(7));
   ASSERT_EQ(-1, fsgrid_detail::linearToX(8));
   ASSERT_EQ(0, fsgrid_detail::linearToX(9));
   ASSERT_EQ(0, fsgrid_detail::linearToX(10));
   ASSERT_EQ(0, fsgrid_detail::linearToX(11));
   ASSERT_EQ(0, fsgrid_detail::linearToX(12));
   ASSERT_EQ(0, fsgrid_detail::linearToX(13));
   ASSERT_EQ(0, fsgrid_detail::linearToX(14));
   ASSERT_EQ(0, fsgrid_detail::linearToX(15));
   ASSERT_EQ(0, fsgrid_detail::linearToX(16));
   ASSERT_EQ(0, fsgrid_detail::linearToX(17));
   ASSERT_EQ(1, fsgrid_detail::linearToX(18));
   ASSERT_EQ(1, fsgrid_detail::linearToX(19));
   ASSERT_EQ(1, fsgrid_detail::linearToX(20));
   ASSERT_EQ(1, fsgrid_detail::linearToX(21));
   ASSERT_EQ(1, fsgrid_detail::linearToX(22));
   ASSERT_EQ(1, fsgrid_detail::linearToX(23));
   ASSERT_EQ(1, fsgrid_detail::linearToX(24));
   ASSERT_EQ(1, fsgrid_detail::linearToX(25));
   ASSERT_EQ(1, fsgrid_detail::linearToX(26));
}

TEST(FsGridTest, linearToY) {
   ASSERT_EQ(-1, fsgrid_detail::linearToY(0));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(1));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(2));
   ASSERT_EQ(0, fsgrid_detail::linearToY(3));
   ASSERT_EQ(0, fsgrid_detail::linearToY(4));
   ASSERT_EQ(0, fsgrid_detail::linearToY(5));
   ASSERT_EQ(1, fsgrid_detail::linearToY(6));
   ASSERT_EQ(1, fsgrid_detail::linearToY(7));
   ASSERT_EQ(1, fsgrid_detail::linearToY(8));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(9));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(10));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(11));
   ASSERT_EQ(0, fsgrid_detail::linearToY(12));
   ASSERT_EQ(0, fsgrid_detail::linearToY(13));
   ASSERT_EQ(0, fsgrid_detail::linearToY(14));
   ASSERT_EQ(1, fsgrid_detail::linearToY(15));
   ASSERT_EQ(1, fsgrid_detail::linearToY(16));
   ASSERT_EQ(1, fsgrid_detail::linearToY(17));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(18));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(19));
   ASSERT_EQ(-1, fsgrid_detail::linearToY(20));
   ASSERT_EQ(0, fsgrid_detail::linearToY(21));
   ASSERT_EQ(0, fsgrid_detail::linearToY(22));
   ASSERT_EQ(0, fsgrid_detail::linearToY(23));
   ASSERT_EQ(1, fsgrid_detail::linearToY(24));
   ASSERT_EQ(1, fsgrid_detail::linearToY(25));
   ASSERT_EQ(1, fsgrid_detail::linearToY(26));
}

TEST(FsGridTest, linearToZ) {
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(0));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(1));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(2));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(3));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(4));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(5));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(6));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(7));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(8));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(9));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(10));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(11));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(12));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(13));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(14));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(15));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(16));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(17));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(18));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(19));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(20));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(21));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(22));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(23));
   ASSERT_EQ(-1, fsgrid_detail::linearToZ(24));
   ASSERT_EQ(0, fsgrid_detail::linearToZ(25));
   ASSERT_EQ(1, fsgrid_detail::linearToZ(26));
}

TEST(FsGridTest, xyzToLinearToxyz) {
   for (uint32_t i = 0; i < 27; i++) {
      const auto x = fsgrid_detail::linearToX(i);
      const auto y = fsgrid_detail::linearToY(i);
      const auto z = fsgrid_detail::linearToZ(i);
      ASSERT_EQ(i, fsgrid_detail::xyzToLinear(x, y, z));
   }
}

TEST(FsGridTest, allNeighboursAreSelf) {
   constexpr auto rank = 1;
   // 13 is never self
   constexpr std::array<int32_t, 27> ranks = {
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };

   for (auto i = 0u; i < 13u; i++) {
      ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[i], 1);
   }

   ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[13], 0);

   for (auto i = 14u; i < 27u; i++) {
      ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[i], 1);
   }
}

TEST(FsGridTest, noneAreSelf) {
   constexpr auto rank = 0;
   // 13 is never self
   constexpr std::array<int32_t, 27> ranks = {
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   };

   for (auto i = 0u; i < 27u; i++) {
      ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[i], 0);
   }
}

TEST(FsGridTest, firstAndLastAreSelf) {
   constexpr auto rank = 0;
   constexpr std::array<int32_t, 27> ranks = {
       0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
   };

   ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[0], 1);
   for (auto i = 1u; i < 26u; i++) {
      ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[i], 0);
   }
   ASSERT_EQ(fsgrid_detail::makeNeigbourBitMask(rank, ranks)[26], 1);
}
