#include <fstream>
#include <iostream>
#include "simulator.hpp"

#define FLOAT            float
#define DOUBLE           double
#define FAST_FIXED(N, K) types::FastFixed<N, K>
#define FIXED(N, K)      types::Fixed<N, K>

#ifndef TYPES
#define TYPES
#endif

// using triple = std::tuple<std::string, std::string, std::string>;
// std::unordered_map<triple, std::function<void()> > registry;

constexpr size_t P = 32, Q = 16;

template<typename FixedCls>
void load_constants(int &n, int &m, int &k, FixedCls rho[256], FixedCls &g, FieldStorageType &field) {
  double g_float;

  std::ifstream input_file("input.txt");

  if (!input_file.is_open()) {
    std::cerr << "Failed to open input file" << std::endl;
    exit(1);
  }
  input_file >> n >> m >> g_float >> k;
  g = FixedCls(g_float);
  input_file.ignore(std::numeric_limits<size_t>::max(), '\n');
  input_file.get();

  field.resize(n);

  for (int i = 0; i < k; i++) {
    char c;
    double f;
    input_file.get(c);
    input_file.get();
    input_file >> f;
    rho[static_cast<int>(c)] = FixedCls(f);
    input_file.ignore(std::numeric_limits<size_t>::max(), '\n');
    input_file.get();
  }

  for (int i = 0; i < n; i++) {
    std::string row;
    std::getline(input_file, row, '\n');
    field[i] = row;
  }
}

template<typename PType, typename VType, typename VFlowType>
void process_type() {
  int n, m;
  int k;
  PType rho[256];
  PType g;
  FieldStorageType field;

  load_constants(n, m, k, rho, g, field);
  Simulator<PType, VType, VFlowType> simulator(n, m, g, rho, field);
  simulator.execute();
}

int main(int argc, char **argv) {
  std::unordered_map<std::string, std::string> arg_map;
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " --p-type=... --v-type=... --v-flow-type=...\n";
    return 1;
  }
  for (int i = 1; i < argc; i++) {
    std::string arg_str = argv[i];
    size_t eq_pos = arg_str.find('=');
    std::string key = arg_str.substr(0, eq_pos);
    std::string value = arg_str.substr(eq_pos + 1);
    arg_map[key] = value;
  }

  std::string p_type, v_type, vf_type;

  if (arg_map.find("--p-type") == arg_map.end() || arg_map.find("--v-type") == arg_map.end() || arg_map.
    find("--v-flow-type") == arg_map.end()) {
    std::cerr << "Missing required argument\n";
    return 1;
  }
  // p_type = arg_map["--p-type"];
  // v_type = arg_map["--v-type"];
  // vf_type = arg_map["--v-flow-type"];

  process_type<Fixed<P, Q, true>, Fixed<P, Q, true>, Fixed<P, Q, true> >();
  return 0;
}
