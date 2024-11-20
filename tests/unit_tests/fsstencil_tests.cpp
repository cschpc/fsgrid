#include <gtest/gtest.h>

#include "fsstencil.hpp"

TEST(FsStencilTest, belowLocal) {
   const fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 0, 0, 0});
   const auto dm = stencil.localityMultipliers({-1, -1, -1});
   ASSERT_EQ(dm[0], -1);
   ASSERT_EQ(dm[1], -1);
   ASSERT_EQ(dm[2], -1);
}

TEST(FsStencilTest, aboveLocal) {
   const fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 0, 0, 0});
   const auto dm = stencil.localityMultipliers({10, 10, 10});
   ASSERT_EQ(dm[0], 1);
   ASSERT_EQ(dm[1], 1);
   ASSERT_EQ(dm[2], 1);
}

TEST(FsStencilTest, insideLocal) {
   const fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 0, 0, 0});
   const auto dm = stencil.localityMultipliers({0, 0, 0});
   ASSERT_EQ(dm[0], 0);
   ASSERT_EQ(dm[1], 0);
   ASSERT_EQ(dm[2], 0);
}

TEST(FsStencilTest, neighbourIndex) {
   const fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 0, 0, 0});
   constexpr std::array xs = {-1, 0, 1};
   constexpr std::array ys = {-1, 0, 1};
   constexpr std::array zs = {-1, 0, 1};

   auto i = 0;
   for (const auto x : xs) {
      for (const auto y : ys) {
         for (const auto z : zs) {
            ASSERT_EQ(stencil.neighbourIndex({x, y, z}), i++);
         }
      }
   }
}

TEST(FsStencilTest, shiftOffsets) {
   constexpr fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 0, 0u, 0});

   constexpr std::array xs = {-1, 0, 1};
   constexpr std::array ys = {-1, 0, 1};
   constexpr std::array zs = {-1, 0, 1};

   constexpr int32_t values[27][3] = {
       {10, 10, 10}, {10, 10, 0},    {10, 10, -10}, {10, 0, 10},    {10, 0, 0},    {10, 0, -10},    {10, -10, 10},
       {10, -10, 0}, {10, -10, -10}, {0, 10, 10},   {0, 10, 0},     {0, 10, -10},  {0, 0, 10},      {0, 0, 0},
       {0, 0, -10},  {0, -10, 10},   {0, -10, 0},   {0, -10, -10},  {-10, 10, 10}, {-10, 10, 0},    {-10, 10, -10},
       {-10, 0, 10}, {-10, 0, 0},    {-10, 0, -10}, {-10, -10, 10}, {-10, -10, 0}, {-10, -10, -10},
   };

   auto i = 0;
   for (const auto x : xs) {
      for (const auto y : ys) {
         for (const auto z : zs) {
            const auto so = stencil.shiftOffsets({x, y, z});
            ASSERT_EQ(so[0], values[i][0]);
            ASSERT_EQ(so[1], values[i][1]);
            ASSERT_EQ(so[2], values[i++][2]);
         }
      }
   }
}

TEST(FsStencilTest, applyMultipliersAndOffset1) {
   constexpr fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {1, 2, 3}, 10, 0u, 0u});
   constexpr std::array value = {6, 7, 8};
   ASSERT_EQ(stencil.applyMultipliersAndOffset(value), 54ul);
}

TEST(FsStencilTest, applyMultipliersAndOffset2) {
   constexpr fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {0, 0, 0}, 10, 0u, 0u});
   constexpr std::array value = {6, 7, 8};
   ASSERT_EQ(stencil.applyMultipliersAndOffset(value), 10ul);
}

TEST(FsStencilTest, applyMultipliersAndOffset3) {
   constexpr fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {1, 1, 1}, 0, 0u, 0u});
   constexpr std::array value = {6, 7, 8};
   ASSERT_EQ(stencil.applyMultipliersAndOffset(value), 21ul);
}

TEST(FsStencilTest, applyMultipliersAndOffset4) {
   constexpr fsgrid::FsStencil stencil = fsgrid::FsStencil({0, 0, 0}, {{10, 10, 10}, {1, 0, 0}, 0, 0u, 0u});
   constexpr std::array value = {6, 7, 8};
   ASSERT_EQ(stencil.applyMultipliersAndOffset(value), 6ul);
}
