#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include "fixed.hpp"

using namespace std;

constexpr size_t T = 1'000'000;

template<typename PType, typename VType, typename VFlowType>
class Simulator {
  public:
    explicit Simulator(
      int n,
      int m,
      PType g,
      const PType _rho[256],
      FieldStorageType &_field
    ): n(n), m(m), g(g) {
      p = std::vector<std::vector<PType> >(n, std::vector<PType>(m));
      old_p = std::vector<std::vector<PType> >(n, std::vector<PType>(m));
      last_use = std::vector<std::vector<int> >(n, std::vector<int>(m));
      dirs = std::vector<std::vector<int> >(n, std::vector<int>(m));
      field = std::move(_field);
      std::memcpy(rho, _rho, sizeof rho);

      velocity = VectorField<VType>(n, m);
      velocity_flow = VectorField<VFlowType>(n, m);

      for (int i = 0; i < n; i++) {
        std::string row;
        for (int j = 0; j < m; j++) {
          char c = field[i][j];
          row.push_back(c);
        }
        this->field.push_back(row);
      }

      rnd.seed(1337);
    }

    tuple<PType, bool, pair<int, int> > propagate_flow(int x, int y, PType lim) {
      last_use[x][y] = UT - 1;
      PType ret = PType{0};
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
          auto cap = velocity.get(x, y, dx, dy);
          auto flow = velocity_flow.get(x, y, dx, dy);
          if (flow == cap) {
            continue;
          }
          auto vp = min(lim, cap - flow);
          if (last_use[nx][ny] == UT - 1) {
            velocity_flow.add(x, y, dx, dy, vp);
            last_use[x][y] = UT;
            return {vp, 1, {nx, ny}};
          }
          auto [t, prop, end] = propagate_flow(nx, ny, vp);
          ret += t;
          if (prop) {
            velocity_flow.add(x, y, dx, dy, t);
            last_use[x][y] = UT;
            return {t, end != pair(x, y), end};
          }
        }
      }
      last_use[x][y] = UT;
      return {ret, 0, {0, 0}};
    }

    void propagate_stop(int x, int y, bool force = false) {
      if (!force) {
        bool stop = true;
        for (auto [dx, dy] : deltas) {
          int nx = x + dx, ny = y + dy;
          if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > VType(0)) {
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
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > VType(0)) {
          continue;
        }
        propagate_stop(nx, ny);
      }
    }

    VType move_prob(int x, int y) {
      VType sum = VType(0);
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < VType(0)) {
          continue;
        }
        sum += v;
      }
      return sum;
    }

    struct ParticleParams {
      char type;
      PType cur_p;
      array<VType, deltas.size()> v;
    };

    void swap_with(ParticleParams &pp, int x, int y) {
      swap(field[x][y], pp.type);
      swap(p[x][y], pp.cur_p);
      swap(velocity.v[x][y], pp.v);
    }

    bool propagate_move(int x, int y, bool is_first) {
      last_use[x][y] = UT - is_first;
      bool ret = false;
      int nx = -1, ny = -1;
      do {
        std::array<VType, deltas.size()> tres;
        VType sum(0);
        for (size_t i = 0; i < deltas.size(); ++i) {
          auto [dx, dy] = deltas[i];
          int nx = x + dx, ny = y + dy;
          if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
            tres[i] = sum;
            continue;
          }
          auto v = velocity.get(x, y, dx, dy);
          if (v < VType(0)) {
            tres[i] = sum;
            continue;
          }
          sum += v;
          tres[i] = sum;
        }

        if (sum == VType(0)) {
          break;
        }

        auto p = PType(rnd) * sum;
        size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

        auto [dx, dy] = deltas[d];
        nx = x + dx;
        ny = y + dy;
        assert(velocity.get(x, y, dx, dy) > VType(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

        ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
      } while (!ret);
      last_use[x][y] = UT;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < VType(0)) {
          propagate_stop(nx, ny);
        }
      }
      if (ret) {
        if (!is_first) {
          ParticleParams pp{};
          swap_with(pp, x, y);
          swap_with(pp, nx, ny);
          swap_with(pp, x, y);
        }
      }
      return ret;
    }

    void execute() {
      // Dir nums for each cell
      for (size_t x = 0; x < n; ++x) {
        for (size_t y = 0; y < m; ++y) {
          if (field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            dirs[x][y] += (field[x + dx][y + dy] != '#');
          }
        }
      }

      // each tick
      for (size_t i = 0; i < T; ++i) {
        PType total_delta_p{};

        // add gravitational force to each velocity
        for (size_t x = 0; x < n; ++x) {
          for (size_t y = 0; y < m; ++y) {
            if (field[x][y] == '#')
              continue;
            if (field[x + 1][y] != '#')
              velocity.add(x, y, 1, 0, g);
          }
        }

        old_p = p;
        for (size_t x = 0; x < n; ++x) {
          for (size_t y = 0; y < m; ++y) {
            if (field[x][y] == '#')
              continue;
            for (auto [dx, dy] : deltas) {
              // Add forces from p
              int nx = x + dx, ny = y + dy;
              if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                auto force = old_p[x][y] - old_p[nx][ny];
                auto &contr = velocity.get(nx, ny, -dx, -dy);
                if (contr * rho[(int) field[nx][ny]] >= force) {
                  contr -= force / rho[(int) field[nx][ny]];
                  continue;
                }
                force -= contr * rho[(int) field[nx][ny]];
                velocity.add(x, y, dx, dy, force / rho[field[x][y]]);
                p[x][y] -= force / PType(dirs[x][y]);
                total_delta_p -= force / PType(dirs[x][y]);
              }
            }
          }
        }

        // Propagate flow
        velocity_flow = VectorField<VFlowType>(n, m);
        bool prop = false;
        do {
          UT += 2;
          prop = false;
          for (size_t x = 0; x < n; ++x) {
            for (size_t y = 0; y < m; ++y) {
              if (field[x][y] != '#' && last_use[x][y] != UT) {
                auto [t, local_prop, _] = propagate_flow(x, y, PType(1));
                if (t > PType(0)) {
                  prop = true;
                }
              }
            }
          }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < n; ++x) {
          for (size_t y = 0; y < m; ++y) {
            if (field[x][y] == '#')
              continue;
            for (auto [dx, dy] : deltas) {
              auto old_v = velocity.get(x, y, dx, dy);
              auto new_v = velocity_flow.get(x, y, dx, dy);
              if (old_v > VType(0)) {
                assert(new_v <= old_v);
                velocity.get(x, y, dx, dy) = new_v;
                auto force = (old_v - new_v) * rho[(int) field[x][y]];
                if (field[x][y] == '.')
                  force *= PType(0.8);
                if (field[x + dx][y + dy] == '#') {
                  p[x][y] += force / PType(dirs[x][y]);
                  total_delta_p += force / PType(dirs[x][y]);
                } else {
                  p[x + dx][y + dy] += force / PType(dirs[x + dx][y + dy]);
                  total_delta_p += force / PType(dirs[x + dx][y + dy]);
                }
              }
            }
          }
        }

        UT += 2;
        prop = false;
        for (size_t x = 0; x < n; ++x) {
          for (size_t y = 0; y < m; ++y) {
            if (field[x][y] != '#' && last_use[x][y] != UT) {
              if (VType(rnd) < VType(move_prob(x, y))) {
                prop = true;
                propagate_move(x, y, true);
              } else {
                propagate_stop(x, y, true);
              }
            }
          }
        }

        if (prop) {
          cout << "Tick " << i << ":\n";
          for (size_t x = 0; x < n; ++x) {
            cout << field[x] << "\n";
          }
        }
      }
    }

  private:
    int n, m;

    FieldStorageType field;
    VectorField<VType> velocity;
    VectorField<VFlowType> velocity_flow;

    PType rho[256];

    std::vector<std::vector<PType> > p, old_p;
    std::vector<std::vector<int> > last_use, dirs;

    int UT = 0;
    PType g = 0.1;

    std::mt19937 rnd;
};
#endif // SIMULATOR_HPP
