#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include <utility>
using namespace std;

constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000;
constexpr array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

char field[N][M + 1] = {
  "####################################################################################",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                       .........                                  #",
  "#..............#            #           .........                                  #",
  "#..............#            #           .........                                  #",
  "#..............#            #           .........                                  #",
  "#..............#            #                                                      #",
  "#..............#            #                                                      #",
  "#..............#            #                                                      #",
  "#..............#            #                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............#                                                      #",
  "#..............#............################                     #                 #",
  "#...........................#....................................#                 #",
  "#...........................#....................................#                 #",
  "#...........................#....................................#                 #",
  "##################################################################                 #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "#                                                                                  #",
  "####################################################################################",
};

template<int N, int K>
struct Fixed {
  static_assert(N > K, "N must be greater than K");
  static_assert(N <= 64, "Maximum supported size is 64 bits");

  using StorageType = typename std::conditional_t<
    (N <= 8), int8_t,
    typename std::conditional_t<
      (N <= 16), int16_t,
      typename std::conditional_t<
        (N <= 32), int32_t,
        int64_t
      >
    >
  >;

  constexpr Fixed(int v) : v(static_cast<StorageType>(v) << K) {
  }
  constexpr Fixed(float f) : v(static_cast<StorageType>(std::round(f * (1 << K)))) {
  }
  constexpr Fixed(double f) : v(static_cast<StorageType>(std::round(f * (1 << K)))) {
  }
  constexpr Fixed() : v(0) {
  }

  static constexpr Fixed from_raw(StorageType x) {
    Fixed ret;
    ret.v = x;
    return ret;
  }

  constexpr double to_double() const {
    return static_cast<double>(v) / (1 << K);
  }

  constexpr float to_float() const {
    return static_cast<float>(v) / (1 << K);
  }

  StorageType v;

  auto operator<=>(const Fixed &) const = default;
  bool operator==(const Fixed &) const = default;

  // Arithmetic operators
  constexpr Fixed operator+(const Fixed &other) const {
    return Fixed::from_raw(v + other.v);
  }

  constexpr Fixed operator-(const Fixed &other) const {
    return Fixed::from_raw(v - other.v);
  }

  constexpr Fixed operator*(const Fixed &other) const {
    return Fixed::from_raw((v * other.v) >> K);
  }

  constexpr Fixed operator/(const Fixed &other) const {
    return Fixed::from_raw((v << K) / other.v);
  }

  Fixed operator+(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v + b.v);
  }

  Fixed operator-(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v - b.v);
  }

  Fixed operator*(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v * b.v) >> 16);
  }

  Fixed operator/(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v << 16) / b.v);
  }

  Fixed &operator+=(Fixed &a, Fixed b) {
    return a = a + b;
  }

  Fixed &operator-=(Fixed &a, Fixed b) {
    return a = a - b;
  }

  Fixed &operator*=(Fixed &a, Fixed b) {
    return a = a * b;
  }

  Fixed &operator/=(Fixed &a, Fixed b) {
    return a = a / b;
  }

  Fixed operator-(Fixed x) {
    return Fixed::from_raw(-x.v);
  }

  Fixed abs(Fixed x) {
    if (x.v < 0) {
      x.v = -x.v;
    }
    return x;
  }

  ostream &operator<<(ostream &out, Fixed x) {
    return out << x.v / (double) (1 << 16);
  }

};



template<int N, int K>
class Simulator {
  using Fixed = Fixed<N, K>;

  static constexpr Fixed inf = Fixed::from_raw(std::numeric_limits<int32_t>::max());
  static constexpr Fixed eps = Fixed::from_raw(deltas.size());

  Fixed operator+(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v + b.v);
  }

  Fixed operator-(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v - b.v);
  }

  Fixed operator*(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v * b.v) >> 16);
  }

  Fixed operator/(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v << 16) / b.v);
  }

  Fixed &operator+=(Fixed &a, Fixed b) {
    return a = a + b;
  }

  Fixed &operator-=(Fixed &a, Fixed b) {
    return a = a - b;
  }

  Fixed &operator*=(Fixed &a, Fixed b) {
    return a = a * b;
  }

  Fixed &operator/=(Fixed &a, Fixed b) {
    return a = a / b;
  }

  Fixed operator-(Fixed x) {
    return Fixed::from_raw(-x.v);
  }

  Fixed abs(Fixed x) {
    if (x.v < 0) {
      x.v = -x.v;
    }
    return x;
  }

  ostream &operator<<(ostream &out, Fixed x) {
    return out << x.v / (double) (1 << 16);
  }

  Fixed rho[256];

  Fixed p[N][M]{}, old_p[N][M];

  struct VectorField {
    array<Fixed, deltas.size()> v[N][M];
    Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
      return get(x, y, dx, dy) += dv;
    }

    Fixed &get(int x, int y, int dx, int dy) {
      size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
      assert(i < deltas.size());
      return v[x][y][i];
    }
  };

  VectorField velocity{}, velocity_flow{};
  int last_use[N][M]{};
  int UT = 0;

  mt19937 rnd(1337);

  tuple<Fixed, bool, pair<int, int> > propagate_flow(int x, int y, Fixed lim) {
    last_use[x][y] = UT - 1;
    Fixed ret = 0;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
        auto cap = velocity.get(x, y, dx, dy);
        auto flow = velocity_flow.get(x, y, dx, dy);
        if (flow == cap) {
          continue;
        }
        // assert(v >= velocity_flow.get(x, y, dx, dy));
        auto vp = min(lim, cap - flow);
        if (last_use[nx][ny] == UT - 1) {
          velocity_flow.add(x, y, dx, dy, vp);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
          return {vp, 1, {nx, ny}};
        }
        auto [t, prop, end] = propagate_flow(nx, ny, vp);
        ret += t;
        if (prop) {
          velocity_flow.add(x, y, dx, dy, t);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
          return {t, prop && end != pair(x, y), end};
        }
      }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
  }

  Fixed random01() {
    return Fixed::from_raw((rnd() & ((1 << 16) - 1)));
  }

  void propagate_stop(int x, int y, bool force = false) {
    if (!force) {
      bool stop = true;
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
          stop = false;
          break;
        }
      }
      if (!stop) {
        return;
      }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
        continue;
      }
      propagate_stop(nx, ny);
    }
  }

  Fixed move_prob(int x, int y) {
    Fixed sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
        continue;
      }
      auto v = velocity.get(x, y, dx, dy);
      if (v < 0) {
        continue;
      }
      sum += v;
    }
    return sum;
  }

  struct ParticleParams {
    char type;
    Fixed cur_p;
    array<Fixed, deltas.size()> v;

    void swap_with(int x, int y) {
      swap(field[x][y], type);
      swap(p[x][y], cur_p);
      swap(velocity.v[x][y], v);
    }
  };

  bool propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
      std::array<Fixed, deltas.size()> tres;
      Fixed sum = 0;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          tres[i] = sum;
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < 0) {
          tres[i] = sum;
          continue;
        }
        sum += v;
        tres[i] = sum;
      }

      if (sum == 0) {
        break;
      }

      Fixed p = random01() * sum;
      size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

      auto [dx, dy] = deltas[d];
      nx = x + dx;
      ny = y + dy;
      assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

      ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
        propagate_stop(nx, ny);
      }
    }
    if (ret) {
      if (!is_first) {
        ParticleParams pp{};
        pp.swap_with(x, y);
        pp.swap_with(nx, ny);
        pp.swap_with(x, y);
      }
    }
    return ret;
  }

  int dirs[N][M]{};
};

