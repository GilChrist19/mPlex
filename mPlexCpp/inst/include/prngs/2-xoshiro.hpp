/* --------------------------------------------------------------------------------
#
#   MGDrivE: declaration of xoshiro class of prngs
#   Marshall Lab
#   Sean Wu (slwu89@berkeley.edu)
#   December 2019
#
-------------------------------------------------------------------------------- */

#ifndef XOSHIRO256ss
#define XOSHIRO256ss


#include <random>
#include <array>
#include <iostream>
#include <limits>
#include <random>


/* --------------------------------------------------------------------------------
#   xoshiro**
-------------------------------------------------------------------------------- */

// for details on what a URNG compatible with STL must have
// see: https://www.boost.org/doc/libs/1_72_0/doc/html/boost_random/reference.html
constexpr std::uint64_t rotl(std::uint64_t x, int k) noexcept
{
  return (x << k) | (x >> (64 - k));
}

class xoshiro256ss
{
public:
  using result_type = std::uint64_t;


  /// \return the smallest value that `operator()` may return. The value is
  ///         strictly less than `max()`
  constexpr static result_type min() noexcept { return 0; }

  /// \return the largest value that `operator()` may return. The value is
  //          strictly greater than `min()`
  constexpr static result_type max() noexcept
  { return std::numeric_limits<result_type>::max(); }

  /// \return a value in the closed interval `[min(), max()]`. Has amortized
  ///         constant complexity
  result_type operator()() noexcept
  {
    const auto result_starstar(rotl(state[1] * 5, 7) * 9);
    const auto t(state[1] << 17);

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= t;

    state[3] = rotl(state[3], 45);

    return result_starstar;
  }

  // set seed
  void seed(const std::array<std::uint64_t, 4> &s) noexcept {
    state = s;
  };

  bool operator==(const xoshiro256ss &rhs) const noexcept
  { return state == rhs.state; }
  bool operator!=(const xoshiro256ss &rhs) const noexcept
  { return !(*this == rhs); }

  friend std::ostream &operator<<(std::ostream &, const xoshiro256ss &);
  friend std::istream &operator>>(std::istream &, xoshiro256ss &);

private:
  std::array<std::uint64_t, 4> state;
};

#endif
