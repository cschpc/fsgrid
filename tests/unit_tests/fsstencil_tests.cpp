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

TEST(FsStencilTest, cellExistsWhenFallBackBitsAreZero) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 0, 0);
   constexpr fsgrid::FsStencil s(0, 0, 0, sc);

   for (int32_t x = -1; x < 2; x++) {
      for (int32_t y = -1; y < 2; y++) {
         for (int32_t z = -1; z < 2; z++) {
            ASSERT_TRUE(s.cellExists(x, y, z));
         }
      }
   }
}

TEST(FsStencilTest, onlyCenterExistsWhenAllFallbackBitsButCenterAreOne) {
   constexpr fsgrid::StencilConstants sc({1, 1, 1}, {0, 0, 0}, 0, 0, 0b00000111111111111101111111111111);
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
