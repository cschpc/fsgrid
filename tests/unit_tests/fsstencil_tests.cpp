#include <gtest/gtest.h>

#include "fsstencil.hpp"

TEST(BitMaskTest, unsetMask) {
   constexpr fsgrid::BitMask32 mask(0);
   for (uint32_t i = 0; i < 32; i++) {
      ASSERT_EQ(mask[i], 0);
   }
}

TEST(BitMaskTest, bit1IsSet) {
   constexpr fsgrid::BitMask32 mask(1);
   ASSERT_EQ(mask[0], 1);
   for (uint32_t i = 1; i < 32; i++) {
      ASSERT_EQ(mask[i], 0);
   }
}

TEST(BitMaskTest, bits1and2AreSet) {
   constexpr fsgrid::BitMask32 mask(3);
   ASSERT_EQ(mask[0], 1);
   ASSERT_EQ(mask[1], 1);
   for (uint32_t i = 2; i < 32; i++) {
      ASSERT_EQ(mask[i], 0);
   }
}

TEST(BitMaskTest, allBitsAreSet) {
   constexpr fsgrid::BitMask32 mask(~0u);
   for (uint32_t i = 0; i < 32; i++) {
      ASSERT_EQ(mask[i], 1);
   }
}

TEST(BitMaskTest, tooLargeIndexGivesZero) {
   constexpr fsgrid::BitMask32 mask(~0u);
   ASSERT_EQ(mask[32], 0);
}

TEST(FsStencilTest, cellExistsWhenFallBackBitsAreZeroAndNumGhostCellsIsOne) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 1, 0, 0);
   constexpr fsgrid::FsStencil s(0, 0, 0, sc);

   for (int32_t x = -1; x < 2; x++) {
      for (int32_t y = -1; y < 2; y++) {
         for (int32_t z = -1; z < 2; z++) {
            ASSERT_TRUE(s.cellExists(x, y, z));
         }
      }
   }
}

TEST(FsStencilTest, cellDoesNotExistWhenFallBackBitsAreZeroAndNumGhostCellsIsZero) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 0, 0, 0);
   constexpr fsgrid::FsStencil s(0, 0, 0, sc);

   for (int32_t x = -1; x < 2; x++) {
      for (int32_t y = -1; y < 2; y++) {
         for (int32_t z = -1; z < 2; z++) {
            if (x == 0 && y == 0 && z == 0) {
               ASSERT_TRUE(s.cellExists(x, y, z));
            } else {
               ASSERT_FALSE(s.cellExists(x, y, z));
            }
         }
      }
   }
}

TEST(FsStencilTest, cellDoesNotExistsWhenFallBackBitsAreZeroAndNumGhostCellsIsOne) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 1, 0, 0);
   constexpr fsgrid::FsStencil s(0, 0, 0, sc);

   for (int32_t x = -2; x < 3; x++) {
      for (int32_t y = -2; y < 3; y++) {
         for (int32_t z = -2; z < 3; z++) {
            if (abs(x) < 2 && abs(y) < 2 && abs(z) < 2) {
               ASSERT_TRUE(s.cellExists(x, y, z));
            } else {
               ASSERT_FALSE(s.cellExists(x, y, z));
            }
         }
      }
   }
}

TEST(FsStencilTest, onlyCenterExistsWhenAllFallbackBitsButCenterAreOne) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 0, 0, 0b00000111111111111101111111111111);
   constexpr fsgrid::FsStencil s(0, 0, 0, sc);

   for (int32_t x = -1; x < 2; x++) {
      for (int32_t y = -1; y < 2; y++) {
         for (int32_t z = -1; z < 2; z++) {
            if (x == 0 && y == 0 && z == 0) {
               ASSERT_TRUE(s.cellExists(x, y, z));
            } else {
               ASSERT_FALSE(s.cellExists(x, y, z));
            }
         }
      }
   }
}

TEST(FsStencilTest, indicesAreCorrect1) {
   // 3x3x3 cube with no ghost cells
   constexpr fsgrid::StencilConstants sc({3, 3, 3}, {1, 3, 9}, 0, 0, 0, 0);
   constexpr fsgrid::FsStencil s(1, 1, 1, sc);

   size_t j = 0;
   for (const auto& i : s.indices()) {
      ASSERT_EQ(i, j++);
   }
}

TEST(FsStencilTest, indicesAreCorrect2) {
   // 3x3x3 cube with 1 ghost cell everywhere, so 5x5x5 cube with ghost cells
   constexpr fsgrid::StencilConstants sc({3, 3, 3}, {1, 5, 25}, 0, 1, 0, 0);
   constexpr fsgrid::FsStencil s(1, 1, 1, sc);

   // clang-format off
   constexpr std::array indices = {
       0, 1, 2,
       5, 6, 7,
       10, 11, 12,
       25, 26, 27,
       30, 31, 32,
       35, 36, 37,
       50, 51, 52,
       55, 56, 57,
       60, 61, 62,
   };
   // clang-format on

   size_t j = 0;
   for (const auto& i : s.indices()) {
      ASSERT_EQ(i, indices[j++]);
   }
}
