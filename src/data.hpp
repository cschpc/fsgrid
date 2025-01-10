#pragma once

/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016-2024 CSC -IT Center for Science

  This file is part of fsgrid

  fsgrid is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  fsgrid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY;
  without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with fsgrid.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <span>

namespace fsgrid {
// Use this as for help, if you need to implement these operations for e.g. Cuda or Hip
struct CMemoryOperations {
   static void* allocate(size_t bytes) { return std::malloc(bytes); }
   static void deallocate(void* ptr) { std::free(ptr); }
   static void memcpy(void* dst, const void* src, size_t bytes, bool = false) { std::memcpy(dst, src, bytes); }
   static void memset(void* dst, int pattern, size_t bytes, bool = false) { std::memset(dst, pattern, bytes); }
};

template <typename T, typename MemOps> struct FsData {
private:
   size_t num_elements = 0ul;
   std::unique_ptr<uint8_t, decltype(&MemOps::deallocate)> memory;

public:
   FsData(size_t num_elements)
       : num_elements(num_elements),
         memory(static_cast<uint8_t*>(MemOps::allocate(getMemReq(num_elements))), &MemOps::deallocate) {
      MemOps::memset(memory.get(), 0, getMemReq(num_elements));
   }

   FsData(const std::span<T> elements) : FsData(elements.size()) {
      MemOps::memcpy(memory.get(), elements.data(), getMemReq(elements.size()));
   }

   [[nodiscard]] static size_t getMemReq(size_t n) { return n * sizeof(T); }

   [[nodiscard]] T& operator[](size_t i) { return data()[i]; }
   [[nodiscard]] const T& operator[](size_t i) const { return data()[i]; }

   [[nodiscard]] size_t size() const { return num_elements; }

   [[nodiscard]] T* data() { return static_cast<T*>(static_cast<void*>(memory.get())); }
   [[nodiscard]] T const* data() const { return static_cast<T*>(static_cast<void*>(memory.get())); }

   [[nodiscard]] std::span<T> view() { return std::span(data(), size()); }
   [[nodiscard]] std::span<const T> view() const { return std::span(data(), size()); }

   void swap(FsData<T, MemOps>& other) noexcept {
      std::swap(num_elements, other.num_elements);
      memory.swap(other.memory);
   }
};

template <typename T, typename MemOps> void swap(FsData<T, MemOps>& a, FsData<T, MemOps>& b) noexcept { a.swap(b); }
} // namespace fsgrid
