#ifndef INCLUDE_MD5SUM
#define INCLUDE_MD5SUM

#include <array>
#include <iterator>
#include <cstdint>
#include "Temp_common.h"
#include "gmpxx.h"

// Code copy pasted from
// https://codereview.stackexchange.com/questions/163872/md5-implementation-in-c11

class md5 {
private:
  std::uint32_t a0_;
  std::uint32_t b0_;
  std::uint32_t c0_;
  std::uint32_t d0_;

  std::array<std::uint32_t, 16> m_array_;
  std::array<std::uint32_t, 16>::iterator m_array_first_;

  static const std::array<std::uint32_t, 64> k_array_;
  static const std::array<std::uint32_t, 64> s_array_;

private:
  static std::uint32_t left_rotate(std::uint32_t x, std::uint32_t c);
       

  template <class OutputIterator>
  static void uint32_to_byte(std::uint32_t n, OutputIterator & first) {

    *first++ = n & 0xff;
    *first++ = (n >> 8) & 0xff;
    *first++ = (n >> 16) & 0xff;
    *first++ = (n >> 24) & 0xff;
  }

  template <class OutputIterator>
  static void uint32_to_hex(std::uint32_t n, OutputIterator & first) {
    const char * hex_chars = "0123456789abcdef";

    std::uint32_t b;

    b = n & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 8) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 16) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 24) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];
  }

private:
  void reset_m_array();

  template <class InputIterator>
  void bytes_to_m_array(InputIterator & first, std::array<std::uint32_t, 16>::iterator m_array_last) {
    for (; m_array_first_ != m_array_last; ++m_array_first_) {
      *m_array_first_ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= *first++ << 16;
      *m_array_first_ |= *first++ << 24;
    }
  }

  template <class InputIterator>
  void true_bit_to_m_array(InputIterator & first, std::ptrdiff_t chunk_length) {
    switch (chunk_length % 4) {
    case 0:
      *m_array_first_++ = 0x00000080;
      break;
    case 1:
      *m_array_first_++ = *first++;
      *m_array_first_ |= 0x00008000;
      break;
    case 2:
      *m_array_first_++ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= 0x00800000;
      break;
    case 3:
      *m_array_first_++ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= *first++ << 16;
      *m_array_first_ |= 0x80000000;
      break;
    }
  }

  void zeros_to_m_array(std::array<std::uint32_t, 16>::iterator m_array_last);

  void original_length_bits_to_m_array(std::uint64_t original_length_bits);

  void hash_chunk();

public:
  template <class InputIterator>
  void update(InputIterator first, InputIterator last) {

    std::uint64_t original_length_bits = std::distance(first, last) * 8;

    std::ptrdiff_t chunk_length;
    while ((chunk_length = std::distance(first, last)) >= 64) {
      reset_m_array();
      bytes_to_m_array(first, m_array_.end());
      hash_chunk();
    }

    reset_m_array();
    bytes_to_m_array(first, m_array_.begin() + chunk_length / 4);
    true_bit_to_m_array(first, chunk_length);

    if (chunk_length >= 56) {
      zeros_to_m_array(m_array_.end());
      hash_chunk();

      reset_m_array();
      zeros_to_m_array(m_array_.end() - 2);
      original_length_bits_to_m_array(original_length_bits);
      hash_chunk();
    }
    else {
      zeros_to_m_array(m_array_.end() - 2);
      original_length_bits_to_m_array(original_length_bits);
      hash_chunk();
    }
  }

public:
  md5();

  template <class Container>
  void digest(Container & container) {
    container.resize(16);
    auto it = container.begin();

    uint32_to_byte(a0_, it);
    uint32_to_byte(b0_, it);
    uint32_to_byte(c0_, it);
    uint32_to_byte(d0_, it);
  }

  template <class Container>
  void hex_digest(Container & container) {
    container.resize(32);
    auto it = container.begin();

    uint32_to_hex(a0_, it);
    uint32_to_hex(b0_, it);
    uint32_to_hex(c0_, it);
    uint32_to_hex(d0_, it);
  }
};

std::string MD5_hash_string(std::string const& data);
mpz_class ConvertHex_to_mpz(std::string const& data);
mpz_class MD5_hash_mpz(std::string const& data);

// Murmurhash function
// Code from wikipedia

static inline uint32_t murmur_32_scramble(uint32_t k) {
  k *= 0xcc9e2d51;
  k = (k << 15) | (k >> 17);
  k *= 0x1b873593;
  return k;
}

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed);

namespace std {
  template <typename T>
  struct hash<std::vector<T>>
  {
    std::size_t operator()(const std::vector<T>& Lval) const
    {
      auto combine_hash=[](size_t & seed, size_t new_hash) -> void {
			  seed ^= new_hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
			};
      int len = Lval.size();
      size_t seed = 0;
      for (int i=0; i<len; i++) {
        size_t e_hash = std::hash<T>()(Lval[i]);
        combine_hash(seed, e_hash);
      }
      return seed;
    }
  };
  template <typename T1, typename T2>
  struct hash<std::pair<T1, T2>>
  {
    std::size_t operator()(const std::pair<T1,T2>& ePair) const
    {
      auto combine_hash=[](size_t & seed, size_t new_hash) -> void {
			  seed ^= new_hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
			};
      size_t seed = std::hash<T1>()(ePair.first);
      size_t e_hash = std::hash<T2>()(ePair.second);
      combine_hash(seed, e_hash);
      return seed;
    }
  };
}




#endif

