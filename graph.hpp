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
  NodeId n;
  EdgeId m;
  sequence<EdgeId> offset;
  sequence<Edge> E;
  Graph() : n(0), m(0) {}
};

struct PchGraph : public Graph {
  EdgeId rm;
  NodeId layer = 0;
  sequence<EdgeId> in_offset;
  sequence<Edge> in_E;
  sequence<NodeId> level;
  sequence<NodeId> rank;
  sequence<NodeId> layerOffset;
  sequence<NodeId> ccOffset;
  sequence<NodeId> ccRank;
};

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
  fprintf(fp, "%u\n", graph.n);
  fprintf(fp, "%zu\n", graph.m);
  fprintf(fp, "%u\n", graph.layer);
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
  for (size_t i = 0; i < graph.layer; i++) {
    fprintf(fp, "%u\n", graph.layerOffset[i]);
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
