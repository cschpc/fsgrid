#include "data.hpp"

#include <gtest/gtest.h>

using namespace fsgrid;

TEST(FsDataTest, size_set_correctly) {
   constexpr size_t N = 10;
   Data<float, CMemoryOperations> data(N);
   ASSERT_EQ(data.size(), N);
}

TEST(FsDataTest, view_works) {
   constexpr size_t N = 10;
   Data<float, CMemoryOperations> data(N);

   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 0.0f);
   }

   for (auto& e : data.view()) {
      e = 666.0f;
   }

   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 666.0f);
   }
}

TEST(FsDataTest, swap_works) {
   constexpr size_t N = 10;
   Data<float, CMemoryOperations> data(N);
   Data<float, CMemoryOperations> data2(2 * N);

   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 0.0f);
   }

   for (const auto& e : data2.view()) {
      ASSERT_EQ(e, 0.0f);
   }

   for (auto& e : data.view()) {
      e = 666.0f;
   }

   using std::swap;
   swap(data, data2);

   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 0.0f);
   }

   for (const auto& e : data2.view()) {
      ASSERT_EQ(e, 666.0f);
   }

   ASSERT_EQ(data.size(), 2 * N);
   ASSERT_EQ(data2.size(), N);
}

TEST(FsDataTest, constructed_correctly_from_elements) {
   constexpr size_t N = 10;
   std::vector<float> elements(N, 1.337f);
   Data<float, CMemoryOperations> data(std::span<float>{elements});
   ASSERT_EQ(data.size(), N);
   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 1.337f);
   }
}

TEST(FsDataTest, data_is_correct) {
   constexpr size_t N = 10;
   Data<size_t, CMemoryOperations> data(N);

   for (const auto& e : data.view()) {
      ASSERT_EQ(e, 0ul);
   }

   size_t i = 0;
   for (auto& e : data.view()) {
      e = i++;
   }

   for (size_t j = 0; j < data.size(); j++) {
      ASSERT_EQ(*(data.data() + j), j);
   }
}

TEST(FsDataTest, indexing_works) {
   constexpr size_t N = 10;
   Data<size_t, CMemoryOperations> data(N);

   size_t i = 0;
   for (auto& e : data.view()) {
      e = i++;
   }

   for (size_t j = 0; j < data.size(); j++) {
      ASSERT_EQ(data.view()[j], data[j]);
   }
}
