#include "hash_functions.h"

std::uint32_t md5::left_rotate(std::uint32_t x, std::uint32_t c) {
  return (x << c) | (x >> (32 - c));
}

void md5::reset_m_array() {
  m_array_first_ = m_array_.begin();
}

void md5::zeros_to_m_array(std::array<std::uint32_t, 16>::iterator m_array_last) {
  for (; m_array_first_ != m_array_last; ++m_array_first_) {
    *m_array_first_ = 0;
  }
}

void md5::original_length_bits_to_m_array(std::uint64_t original_length_bits) {
  original_length_bits &= 0xffffffffffffffff;
  *m_array_first_++ = (original_length_bits) & 0x00000000ffffffff;
  *m_array_first_++ = (original_length_bits & 0xffffffff00000000) >> 32;
}

void md5::hash_chunk() {
  std::uint32_t A = a0_;
  std::uint32_t B = b0_;
  std::uint32_t C = c0_;
  std::uint32_t D = d0_;

  std::uint32_t F;
  unsigned int g;

  for (unsigned int i = 0; i < 64; ++i) {
    if (i < 16) {
      F = (B & C) | ((~B) & D);
      g = i;
    }
    else if (i < 32) {
      F = (D & B) | ((~D) & C);
      g = (5 * i + 1) & 0xf;
    }
    else if (i < 48) {
      F = B ^ C ^ D;
      g = (3 * i + 5) & 0xf;
    }
    else {
      F = C ^ (B | (~D));
      g = (7 * i) & 0xf;
    }

    std::uint32_t D_temp = D;
    D = C;
    C = B;
    B += left_rotate(A + F + k_array_[i] + m_array_[g], s_array_[i]);
    A = D_temp;
  }

  a0_ += A;
  b0_ += B;
  c0_ += C;
  d0_ += D;
}

md5::md5()
  : a0_(0x67452301),
    b0_(0xefcdab89),
    c0_(0x98badcfe),
    d0_(0x10325476)
{}

const std::array<std::uint32_t, 64> md5::k_array_ = {
						     0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee,
						     0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
						     0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
						     0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
						     0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
						     0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
						     0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
						     0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
						     0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
						     0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
						     0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
						     0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
						     0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039,
						     0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
						     0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
						     0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391
};

const std::array<std::uint32_t, 64> md5::s_array_ = {
						     7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,
						     5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,
						     4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,
						     6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21
};


std::string MD5_hash_string(std::string const& data)
{
  std::string data_hex_digest;
  md5 hash;
  hash.update(data.begin(), data.end());
  hash.hex_digest(data_hex_digest);
  //
  std::stringstream s;
  s << data_hex_digest;
  std::string converted(s.str());
  return converted;
}

mpz_class ConvertHex_to_mpz(std::string const& data)
{
  std::vector<std::string> LChar{"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"};
  auto GetPosition=[&](std::string const& eChar) -> int {
		     for (int iPos=0; iPos<16; iPos++)
		       if (LChar[iPos] == eChar)
			 return iPos;
		     std::cerr << "Wrong character\n";
		     throw TerminalException{1};
		   };
  mpz_class sum=0;
  mpz_class epow=1;
  for (size_t u=0; u<data.size(); u++) {
    int pos = GetPosition(data.substr(u,1));
    sum += pos * epow;
    epow *= 16;
  }
  return sum;
}


mpz_class MD5_hash_mpz(std::string const& data)
{
  std::string data_out = MD5_hash_string(data);
  return ConvertHex_to_mpz(data_out);
}


// Murmurhash function
// Code from wikipedia

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed)
{
  uint32_t h = seed;
  uint32_t k;
  /* Read in groups of 4. */
  for (size_t i = len >> 2; i; i--) {
    // Here is a source of differing results across endiannesses.
    // A swap here has no effects on hash properties though.
    memcpy(&k, key, sizeof(uint32_t));
    key += sizeof(uint32_t);
    h ^= murmur_32_scramble(k);
    h = (h << 13) | (h >> 19);
    h = h * 5 + 0xe6546b64;
  }
  /* Read the rest. */
  k = 0;
  for (size_t i = len & 3; i; i--) {
    k <<= 8;
    k |= key[i - 1];
  }
  // A swap is *not* necessary here because the preceding loop already
  // places the low bytes in the low places according to whatever endianness
  // we use. Swaps only apply when the memory is copied in a chunk.
  h ^= murmur_32_scramble(k);
  /* Finalize. */
  h ^= len;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}




