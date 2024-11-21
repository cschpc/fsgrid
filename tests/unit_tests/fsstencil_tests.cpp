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
