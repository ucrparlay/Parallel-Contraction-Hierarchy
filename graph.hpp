#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <fstream>
#include <map>
#include <vector>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/sequence.h"
using namespace std;
using namespace parlay;

#define ELONG
#ifdef NLONG
using NodeId = uint64_t;
#else
using NodeId = uint32_t;
#endif

#ifdef ELONG
using EdgeId = uint64_t;
#else
using EdgeId = uint32_t;
#endif

using EdgeTy = uint32_t;
constexpr NodeId UINT_N_MAX = numeric_limits<NodeId>::max();
constexpr EdgeId UINT_E_MAX = numeric_limits<EdgeId>::max();
constexpr EdgeTy DIST_MAX = numeric_limits<EdgeTy>::max();
constexpr int LOG2_WEIGHT = 5;
constexpr int WEIGHT = 1 << LOG2_WEIGHT;

struct Edge {
  NodeId v;
  EdgeTy w;
  NodeId hop;
  // NodeId mid;
  Edge() : v(0), w(0), hop(1){};
  Edge(NodeId _v, EdgeTy _w) : v(_v), w(_w), hop(1) {}
  bool operator<(const Edge &rhs) const {
    if (v != rhs.v) {
      return v < rhs.v;
    }
    return w < rhs.w;
  }
  bool operator!=(const Edge &rhs) const { return v != rhs.v || w != rhs.w; }
};

struct Graph {
  size_t n;
  size_t m;
  sequence<EdgeId> offset;
  sequence<Edge> E;
  Graph() : n(0), m(0) {}
};

struct PchGraph : public Graph {
  size_t rm;
  size_t layer;
  sequence<EdgeId> in_offset;
  sequence<Edge> in_E;
  sequence<NodeId> level;
  sequence<NodeId> rank;
  sequence<NodeId> layerOffset;
  sequence<NodeId> ccOffset;
  sequence<NodeId> ccRank;

  void read_pch_format(char *const filename);
  void write_pch_format(char *const filename);
};

void PchGraph::write_pch_format(char *const filename) {
  ofstream ofs(filename);
  if (!ofs.is_open()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    abort();
  }
  size_t num_cc = 0;
  ofs.write(reinterpret_cast<char *>(&n), sizeof(size_t));
  ofs.write(reinterpret_cast<char *>(&m), sizeof(size_t));
  ofs.write(reinterpret_cast<char *>(&rm), sizeof(size_t));
  ofs.write(reinterpret_cast<char *>(&layer), sizeof(size_t));
  ofs.write(reinterpret_cast<char *>(&num_cc), sizeof(size_t));
  ofs.write(reinterpret_cast<char *>(offset.begin()), sizeof(EdgeId) * (n + 1));
  ofs.write(reinterpret_cast<char *>(E.begin()), sizeof(Edge) * m);
  ofs.write(reinterpret_cast<char *>(in_offset.begin()),
            sizeof(EdgeId) * (n + 1));
  ofs.write(reinterpret_cast<char *>(in_E.begin()), sizeof(Edge) * m);
  assert(level.size() == n);
  ofs.write(reinterpret_cast<char *>(level.begin()), sizeof(NodeId) * n);
  assert(rank.size() == n);
  ofs.write(reinterpret_cast<char *>(rank.begin()), sizeof(NodeId) * n);
  assert(layerOffset.size() == n + 1);
  ofs.write(reinterpret_cast<char *>(layerOffset.begin()), sizeof(NodeId) * n);
  assert(ccOffset.size() == num_cc + 1);
  ofs.write(reinterpret_cast<char *>(ccOffset.begin()), sizeof(NodeId) * n);
  assert(ccRank.size() == n);
  ofs.write(reinterpret_cast<char *>(ccRank.begin()), sizeof(NodeId) * n);
  ofs.close();
}

void PchGraph::read_pch_format(char *const filename) {
  ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    abort();
  }
  size_t num_cc = 0;
  ifs.read(reinterpret_cast<char *>(&n), sizeof(size_t));
  ifs.read(reinterpret_cast<char *>(&m), sizeof(size_t));
  ifs.read(reinterpret_cast<char *>(&rm), sizeof(size_t));
  ifs.read(reinterpret_cast<char *>(&layer), sizeof(size_t));
  ifs.read(reinterpret_cast<char *>(&num_cc), sizeof(size_t));
  ifs.read(reinterpret_cast<char *>(offset.begin()), sizeof(EdgeId) * (n + 1));
  ifs.read(reinterpret_cast<char *>(E.begin()), sizeof(Edge) * m);
  ifs.read(reinterpret_cast<char *>(in_offset.begin()),
           sizeof(EdgeId) * (n + 1));
  ifs.read(reinterpret_cast<char *>(in_E.begin()), sizeof(Edge) * m);
  assert(level.size() == n);
  ifs.read(reinterpret_cast<char *>(level.begin()), sizeof(NodeId) * n);
  assert(rank.size() == n);
  ifs.read(reinterpret_cast<char *>(rank.begin()), sizeof(NodeId) * n);
  assert(layerOffset.size() == n + 1);
  ifs.read(reinterpret_cast<char *>(layerOffset.begin()), sizeof(NodeId) * n);
  assert(ccOffset.size() == num_cc + 1);
  ifs.read(reinterpret_cast<char *>(ccOffset.begin()), sizeof(NodeId) * n);
  assert(ccRank.size() == n);
  ifs.read(reinterpret_cast<char *>(ccRank.begin()), sizeof(NodeId) * n);
  ifs.close();
}

bool is_space(char c) {
  switch (c) {
    case '\r':
    case '\t':
    case '\n':
    case ' ':
    case 0:
      return true;
    default:
      return false;
  }
}

void generate_weight(Graph &graph) {
  parallel_for(0, graph.n, [&](size_t i) {
    for (size_t j = graph.offset[i]; j < graph.offset[i + 1]; j++) {
      graph.E[j].w =
          ((hash32_2(i) ^ hash32_2(graph.E[j].v)) & (WEIGHT - 1)) + 1;
    }
  });
}

Graph read_binary_format(char const *filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    fprintf(stderr, "Error: Cannot open file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  if (fstat(fd, &sb) == -1) {
    fprintf(stderr, "Error: Unable to acquire file stat\n");
    exit(EXIT_FAILURE);
  }
  char *data =
      static_cast<char *>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  size_t len = sb.st_size;
  size_t n = reinterpret_cast<uint64_t *>(data)[0];
  size_t m = reinterpret_cast<uint64_t *>(data)[1];
  size_t sizes = reinterpret_cast<uint64_t *>(data)[2];
  assert(sizes == (n + 1) * 8 + m * 4 + 3 * 8);
  Graph graph;
  graph.n = n;
  graph.m = m;
  graph.offset = sequence<EdgeId>(n + 1);
  graph.E = sequence<Edge>(m);
  parallel_for(0, n + 1, [&](size_t i) {
    graph.offset[i] = reinterpret_cast<uint64_t *>(data + 3 * 8)[i];
  });
  parallel_for(0, m, [&](size_t i) {
    graph.E[i].v = reinterpret_cast<uint32_t *>(data + 3 * 8 + (n + 1) * 8)[i];
  });
  generate_weight(graph);
  if (data) {
    const void *b = data;
    munmap(const_cast<void *>(b), len);
  }
  return graph;
}

Graph read_pbbs(const char *const filename) {
  ifstream file(filename, ifstream::in | ifstream::binary | ifstream::ate);
  if (!file.is_open()) {
    cerr << "Error: Cannot open file " << filename << '\n';
    abort();
  }
  long end = file.tellg();
  file.seekg(0, file.beg);
  long length = end - file.tellg();

  sequence<char> strings(length + 1);
  file.read(strings.begin(), length);
  file.close();

  sequence<bool> flag(length);
  parallel_for(0, length, [&](size_t i) {
    if (is_space(strings[i])) {
      strings[i] = 0;
    }
  });
  flag[0] = strings[0];
  parallel_for(1, length, [&](size_t i) {
    if (is_space(strings[i - 1]) && !is_space(strings[i])) {
      flag[i] = true;
    } else {
      flag[i] = false;
    }
  });

  auto offset_seq = pack_index(flag);
  sequence<char *> words(offset_seq.size());
  parallel_for(0, offset_seq.size(),
               [&](size_t i) { words[i] = strings.begin() + offset_seq[i]; });

  Graph graph;
  size_t n = atol(words[1]);
  size_t m = atol(words[2]);
  graph.n = n;
  graph.m = m;
  graph.offset = sequence<EdgeId>(n + 1);
  graph.E = sequence<Edge>(m);
  assert(words.size() == n + m + m + 3);
  parallel_for(0, n, [&](size_t i) { graph.offset[i] = atol(words[i + 3]); });
  graph.offset[n] = m;
  parallel_for(0, m, [&](size_t i) { graph.E[i].v = atol(words[n + i + 3]); });
  parallel_for(0, m, [&](size_t i) {
    graph.E[i].w = (EdgeTy)atof(words[n + m + i + 3]);
  });
  return graph;
}

void write_pbbs_format(PchGraph &graph, char const *filename) {
  printf("Info: Writing pbbs format\n");
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "WeightedAdjacency\n");
  fprintf(fp, "%zu\n", graph.n);
  fprintf(fp, "%zu\n", graph.m);
  fprintf(fp, "%zu\n", graph.rm);
  fprintf(fp, "%zu\n", graph.layerOffset.size());
  fprintf(fp, "%zu\n", graph.ccOffset.size());
  for (size_t i = 0; i < graph.n; i++) {
    fprintf(fp, "%lu\n", graph.offset[i]);
  }
  for (size_t i = 0; i < graph.n; i++) {
    fprintf(fp, "%u\n", graph.rank[i]);
  }
  for (size_t i = 0; i < graph.m; i++) {
    fprintf(fp, "%d\n", graph.E[i].v);
  }
  for (size_t i = 0; i < graph.m; i++) {
    fprintf(fp, "%d\n", graph.E[i].w);
  }
  for (size_t i = 0; i < graph.layerOffset.size(); i++) {
    fprintf(fp, "%u\n", graph.layerOffset[i]);
  }
  for (size_t i = 0; i < graph.n; i++) {
    fprintf(fp, "%lu\n", graph.in_offset[i]);
  }
  for (size_t i = 0; i < graph.rm; i++) {
    fprintf(fp, "%d\n", graph.in_E[i].v);
  }
  for (size_t i = 0; i < graph.rm; i++) {
    fprintf(fp, "%d\n", graph.in_E[i].w);
  }
  for (size_t i = 0; i < graph.ccOffset.size(); i++) {
    fprintf(fp, "%u\n", graph.ccOffset[i]);
  }
  for (size_t i = 0; i < graph.n; i++) {
    fprintf(fp, "%d\n", graph.ccRank[i]);
  }
  for (size_t i = 0; i < graph.ccOffset.size() - 1; i++) {
    fprintf(fp, "%d\n",
            (graph.level[graph.layerOffset[graph.ccOffset[i]]] != 0));
  }
  fclose(fp);
}

Graph read_graph(char *filename) {
  string str_filename = string(filename);
  size_t idx = str_filename.find_last_of('.');
  if (idx == string::npos) {
    cerr << "Error: No graph extension provided\n";
    abort();
  }
  string subfix = str_filename.substr(idx + 1);
  Graph graph;
  if (subfix == "adj") {
    graph = read_pbbs(filename);
  } else if (subfix == "bin") {
    graph = read_binary_format(filename);
  } else {
    cerr << "Error: Invalid graph extension\n";
  }
  return graph;
}

Graph generate_random_graph(size_t n, size_t m, int seed = 0) {
  Graph G;
  G.n = n;
  G.m = m;
  constexpr int mul = 1e9 + 7;
  constexpr int MAX_WEIGHT = 1 << 10;
  auto edgelist = parlay::sequence<std::pair<NodeId, NodeId>>::uninitialized(m);
  parlay::parallel_for(0, m, [&](size_t i) {
    NodeId u = parlay::hash32(seed * mul + i * 2) % n;
    NodeId v = parlay::hash32(seed * mul + i * 2 + 1) % n;
    edgelist[i] = std::make_pair(u, v);
  });
  parlay::sort_inplace(parlay::make_slice(edgelist),
                       [](const std::pair<NodeId, NodeId> &a,
                          const std::pair<NodeId, NodeId> &b) {
                         if (a.first != b.first) {
                           return a.first < b.first;
                         }
                         return a.second < b.second;
                       });
  G.offset = parlay::sequence<EdgeId>(n + 1, m);
  G.E = parlay::sequence<Edge>(m);
  parlay::parallel_for(0, m, [&](size_t i) {
    G.E[i].v = edgelist[i].second;
    G.E[i].w = parlay::hash32(seed * mul + i) % MAX_WEIGHT;
    if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
      G.offset[edgelist[i].first] = i;
    }
  });
  parlay::scan_inclusive_inplace(
      parlay::make_slice(G.offset.rbegin(), G.offset.rend()),
      parlay::minm<EdgeId>());
  return G;
}

void read_contracted_graph(const char *const filename, PchGraph &contract) {
  ifstream file(filename, ifstream::in | ifstream::binary | ifstream::ate);
  if (!file.is_open()) {
    cerr << "Error: Cannot open file " << filename << '\n';
    abort();
  }
  long end = file.tellg();
  file.seekg(0, file.beg);
  long length = end - file.tellg();

  sequence<char> strings(length + 1);
  file.read(strings.begin(), length);
  file.close();

  sequence<bool> flag(length);
  parallel_for(0, length, [&](size_t i) {
    if (is_space(strings[i])) {
      strings[i] = 0;
    }
  });
  flag[0] = strings[0];
  parallel_for(1, length, [&](size_t i) {
    if (is_space(strings[i - 1]) && !is_space(strings[i])) {
      flag[i] = true;
    } else {
      flag[i] = false;
    }
  });

  auto offset_seq = pack_index(flag);
  sequence<char *> words(offset_seq.size());
  parallel_for(0, offset_seq.size(),
               [&](size_t i) { words[i] = strings.begin() + offset_seq[i]; });

  size_t n = atol(words[1]);
  contract.n = n;
  contract.m = atol(words[2]);
  contract.rm = atol(words[3]);
  size_t layer = atol(words[4]);
  size_t cc_nums = atol(words[5]);
  assert(words.size() == n + n + n + n + contract.m + contract.m + contract.rm +
                             contract.rm + layer + cc_nums + cc_nums + 5);
  contract.layer = layer;
  contract.offset = sequence<EdgeId>(n + 1);
  contract.in_offset = sequence<EdgeId>(n + 1);
  contract.rank = sequence<NodeId>(n);
  contract.layerOffset = sequence<NodeId>(layer + 1);
  contract.level = sequence<NodeId>(n + 1);
  contract.ccOffset = sequence<NodeId>(cc_nums);
  contract.ccRank = sequence<NodeId>(n);
  contract.E = sequence<Edge>(contract.m);
  contract.in_E = sequence<Edge>(contract.rm);
  sequence<bool> all_contracted(cc_nums);
  parallel_for(0, n,
               [&](size_t i) { contract.offset[i] = atol(words[i + 6]); });
  contract.offset[n] = contract.m;
  parallel_for(0, n,
               [&](size_t i) { contract.rank[i] = atol(words[i + n + 6]); });
  parallel_for(0, contract.m,
               [&](size_t i) { contract.E[i].v = atol(words[i + n + n + 6]); });
  parallel_for(0, contract.m, [&](size_t i) {
    contract.E[i].w = atol(words[i + n + n + contract.m + 6]);
  });
  parallel_for(0, layer, [&](size_t i) {
    contract.layerOffset[i] =
        atol(words[i + n + n + contract.m + contract.m + 6]);
  });
  parallel_for(0, n, [&](size_t i) {
    contract.in_offset[i] =
        atol(words[i + n + n + contract.m + contract.m + layer + 6]);
  });
  contract.in_offset[n] = contract.rm;
  parallel_for(0, contract.rm, [&](size_t i) {
    contract.in_E[i].v =
        atol(words[i + n + n + n + contract.m + contract.m + layer + 6]);
  });
  parallel_for(0, contract.rm, [&](size_t i) {
    contract.in_E[i].w = atol(words[i + n + n + n + contract.m + contract.m +
                                    contract.rm + layer + 6]);
  });
  parallel_for(0, cc_nums, [&](size_t i) {
    contract.ccOffset[i] = atol(words[i + n + n + n + contract.m + contract.m +
                                      contract.rm + contract.rm + layer + 6]);
  });
  parallel_for(0, n, [&](size_t i) {
    contract.ccRank[i] =
        atol(words[i + n + n + n + contract.m + contract.m + contract.rm +
                   contract.rm + layer + cc_nums + 6]);
  });
  parallel_for(0, cc_nums - 1, [&](size_t i) {
    all_contracted[i] =
        atol(words[i + n + n + n + n + contract.m + contract.m + contract.rm +
                   contract.rm + layer + cc_nums + 6]);
  });
  parallel_for(0, cc_nums - 1, [&](size_t i) {
    size_t stp = all_contracted[i];
    parallel_for(contract.ccOffset[i], contract.ccOffset[i + 1], [&](size_t j) {
      parallel_for(contract.layerOffset[j], contract.layerOffset[j + 1],
                   [&](size_t k) {
                     contract.level[k] = stp + j - contract.ccOffset[i];
                   });
    });
  });
}
