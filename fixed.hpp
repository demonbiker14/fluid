//
// Created by Daniel Chiliaev on 21/12/2024.
//

#ifndef TYPES_H
#define TYPES_H

#include <iostream>

constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};



using FieldStorageType = std::vector<std::string>;


template<int P, int Q, bool fast = false>
struct Fixed {
  static_assert(P > Q, "P must be greater than Q");
  static_assert(P <= 64, "Maximum supported size is 64 bits");

  using StorageType = std::conditional_t<
    fast, std::conditional_t<
      (P <= 8), int_fast8_t,
      std::conditional_t<
        (P <= 16), int_fast16_t,
        std::conditional_t<
          (P <= 32), int_fast32_t,
          int_fast64_t
        >
      >
    >, std::conditional_t<
      (P <= 8), int8_t,
      std::conditional_t<
        (P <= 16), int16_t,
        std::conditional_t<
          (P <= 32), int32_t,
          int64_t
        >
      >
    > >;

  constexpr explicit Fixed(std::mt19937 &rnd) : v(rnd() & ((1 << Q) - 1)) {
  }

  explicit constexpr Fixed(int v) : v(static_cast<StorageType>(v) << Q) {
  }
  explicit constexpr Fixed(float f) : v(static_cast<StorageType>(std::round(f * (1 << Q)))) {
  }
  explicit constexpr Fixed(double f) : v(static_cast<StorageType>(std::round(f * (1 << Q)))) {
  }
  constexpr Fixed() : v(0) {
  }

  static constexpr Fixed from_raw(StorageType x) {
    Fixed ret;
    ret.v = x;
    return ret;
  }

  constexpr double to_double() const {
    return static_cast<double>(v) / (1 << Q);
  }

  constexpr float to_float() const {
    return static_cast<float>(v) / (1 << Q);
  }

  auto operator<=>(const Fixed &) const = default;
  bool operator==(const Fixed &) const = default;

  constexpr Fixed operator+(const Fixed &other) const {
    return Fixed::from_raw(v + other.v);
  }

  constexpr Fixed operator-(const Fixed &other) const {
    return Fixed::from_raw(v - other.v);
  }

  constexpr Fixed operator*(const Fixed &other) const {
    return Fixed::from_raw((v * other.v) >> Q);
  }

  constexpr Fixed operator/(const Fixed &other) const {
    return Fixed::from_raw((v << Q) / other.v);
  }

  Fixed &operator=(const int &&num) {
    *this = Fixed::from_raw(num);
    return *this;
  }

  Fixed &operator=(const float &&num) {
    *this = Fixed::from_raw(num);
    return *this;
  }

  Fixed &operator=(const double &&num) {
    *this = Fixed::from_raw(num);
    return *this;
  }

  Fixed &operator+=(const Fixed &other) {
    return *this = *this + other;
  }

  Fixed &operator-=(const Fixed &other) {
    return *this = *this - other;
  }

  Fixed &operator*=(const Fixed &other) {
    return *this = *this * other;
  }

  Fixed &operator/=(const Fixed &other) {
    return *this = *this / other;
  }
  constexpr Fixed operator-() const {
    return Fixed::from_raw(-v);
  }

  constexpr Fixed abs() const {
    return Fixed::from_raw(v < 0 ? -v : v);
  }

  friend std::ostream &operator<<(std::ostream &out, const Fixed &x) {
    return out << x.v / (double) (1 << Q);
  }

  StorageType inf() {
    return std::numeric_limits<StorageType>::max();
  }

  StorageType eps() {
    return from_raw(deltas.size());
  }

  StorageType v;
};

template<typename FixedType>
struct VectorField {
  using Fixed = FixedType;

  std::vector<std::vector<
    std::array<Fixed, deltas.size()>
  > > v;

  VectorField() = default;
  VectorField(int n, int m) {
    v.resize(n, std::vector<std::array<Fixed, deltas.size()> >(m));
  }

  Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
    return get(x, y, dx, dy) += dv;
  }

  Fixed &get(int x, int y, int dx, int dy) {
    size_t i = std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
    assert(i < deltas.size());
    return v[x][y][i];
  }
};


#endif //TYPES_H
