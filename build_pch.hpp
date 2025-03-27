
#include <stdlib.h>

#include <unordered_map>

#include "connectivity.hpp"
#include "graph.hpp"
#include "hash_bag.h"
#include "hash_map.h"
#include "parlay/primitives.h"
#include "utilities.h"

using namespace std;
using namespace parlay;
char *INPUT_FILEPATH = nullptr;
char *OUTPUT_FILEPATH = nullptr;

#define EDGEDIFF_SAMPLES 150
struct Node {
  float edge_diff;
  NodeId removed_neighbor_num;
  NodeId out_degree;
  NodeId in_degree;
  // bool vertices_settled;
  // bool vertices_contracted;
};
struct PCH {
 private:
  bool degree_ordering;
  NodeId g_degree_bound;
  NodeId g_max_pop_count;
  NodeId g_tot_level = 0;
  double g_upper_score_bound = 0;
  double g_sample_bound = 0;
  bool early_stop = false;
  sequence<Node> info;
  sequence<NodeId> overlay_vertices;
  sequence<NodeId> vertices_need_score_update;
  sequence<size_t> priority;
  sequence<bool> vertices_settled;
  sequence<bool> vertices_contracted;

  hash_map<NodeId, NodeId, EdgeTy> forward_edges_in_ch_hash;
  hash_map<NodeId, NodeId, EdgeTy> backward_edges_in_ch_hash;
  hash_map<NodeId, NodeId, EdgeTy> newly_inserted_out_edges_hash;
  hash_map<NodeId, NodeId, EdgeTy> newly_inserted_in_edges_hash;
  hash_map2<NodeId, NodeId, EdgeTy> vertices_around_hash;
  hashbag<uint64_t> bag;
  sequence<uint64_t> neighbors;

  // helper func
  bool check_edge_valid(size_t idx, NodeId v, bool in_csr);
  void transfer_neighbors(NodeId u, vector<NodeId> &idx, vector<NodeId> &hop, vector<EdgeTy> &wgh,
                          bool forward);
  EdgeId clip_from_hash(sequence<Edge> &E,
                        hash_map<NodeId, NodeId, EdgeTy> &hashMap, NodeId nd,
                        EdgeId off, bool judge = true);
  // remove self loop and redundant edges
  void createInEdgesFromOutEdges();
  void preprocess(bool remove_self_loops_and_parallel_edges);
  void initialize();
  void transferHashtoCSR();
  void buildContractionHierarchy();
  void calDisToVerticesAround();
  double sampleUpperBound();
  bool eligibleForRemoval(NodeId s);
  void pruneNeighbors(NodeId s);
  bool calScore();
  bool insertHelper(NodeId left, NodeId right, EdgeTy len, NodeId hop);
  void reorderByLayerAndCC(sequence<pair<NodeId, NodeId>> &node_and_ids);
  void setOrderedLayer();
  void dumpCH();
  void checkDegree();

  template <typename F>
  bool iterateHashTable(NodeId u, F f);

 public:
  const Graph &G_in;
  PchGraph G;
  PchGraph GC;
  bool print_detail;
  PchGraph createContractionHierarchy();
  PCH(Graph &_G_in, int _max_pop_count = 500, bool _degree_ordering = false,
      NodeId _g_degree_bound = 30, double _g_sample_bound = 0.05,
      bool _print_detail = false)
      : degree_ordering(_degree_ordering),
        g_degree_bound(_g_degree_bound),
        g_max_pop_count(_max_pop_count),
        g_sample_bound(_g_sample_bound),
        info(_G_in.n),
        priority(_G_in.n),
        forward_edges_in_ch_hash(2 * _G_in.m),
        backward_edges_in_ch_hash(2 * _G_in.m),
        newly_inserted_out_edges_hash(0.5 * _G_in.m),
        newly_inserted_in_edges_hash(0.5 * _G_in.m),
        vertices_around_hash(20 * _G_in.m),
        bag(5 * _G_in.n),
        neighbors(5 * _G_in.n),
        G_in(_G_in),
        print_detail(_print_detail) {}
};

bool PCH::check_edge_valid([[maybe_unused]] size_t idx, NodeId v, bool in_csr) {
  if (in_csr) {
    if (vertices_contracted[v]) return false;
    return true;
  }
  return false;
}

EdgeId PCH::clip_from_hash(sequence<Edge> &E,
                           hash_map<NodeId, NodeId, EdgeTy> &hashMap, NodeId nd,
                           EdgeId off, bool judge) {
  unsigned long index = hashMap.first_index(nd);
  while (true) {
    if (hashMap.H[index].valid) {
      if (hashMap.H[index].key == UINT_N_MAX) break;
      if (hashMap.H[index].key == nd) {
        if (!judge || !vertices_contracted[hashMap.H[index].value.first]) {
          E[off].v = hashMap.H[index].value.first;
          E[off].w = hashMap.H[index].value.second;
          E[off].hop = hashMap.H[index].hop;
          off++;
        }
      }
    }
    index = hashMap.next_index(index);
  }
  return off;
}

void PCH::initialize() {
  GC.n = G.n;
  GC.m = G.m;
  G.level = sequence<NodeId>(G.n, 1);
  GC.offset = sequence<EdgeId>(GC.n + 1);
  GC.in_offset = sequence<EdgeId>(GC.n + 1);
  GC.E = sequence<Edge>(GC.m);
  GC.in_E = sequence<Edge>(GC.m);
  overlay_vertices = sequence<NodeId>(G.n);
  parallel_for(0, G.n, [&](size_t i) {
    overlay_vertices[i] = i;
    info[i].out_degree = G.offset[i + 1] - G.offset[i];
    info[i].in_degree = G.in_offset[i + 1] - G.in_offset[i];
    priority[i] = i + 1;
    info[i].removed_neighbor_num = 0;
  });
  priority = parlay::random_shuffle(priority);
}

void PCH::createInEdgesFromOutEdges() {
  size_t n = G.n;
  G.rm = G.m;
  sequence<tuple<NodeId, NodeId, EdgeTy>> edgelist(G.m);
  parallel_for(0, n, [&](NodeId u) {
    parallel_for(G.offset[u], G.offset[u + 1], [&](EdgeId i) {
      edgelist[i] = make_tuple(G.E[i].v, u, G.E[i].w);
    });
  });
  sort_inplace(make_slice(edgelist));
  G.in_offset = sequence<EdgeId>(n + 1, G.m);
  G.in_E = sequence<Edge>(G.m);
  parallel_for(0, G.m, [&](size_t i) {
    G.in_E[i].v = get<1>(edgelist[i]);
    G.in_E[i].w = get<2>(edgelist[i]);
    if (i == 0 || get<0>(edgelist[i]) != get<0>(edgelist[i - 1])) {
      G.in_offset[get<0>(edgelist[i])] = i;
    }
  });
  parlay::scan_inclusive_inplace(
      parlay::make_slice(G.in_offset.rbegin(), G.in_offset.rend()),
      parlay::minm<EdgeId>());
}

void PCH::preprocess(bool remove_self_loops_and_parallel_edges = true) {
  if (remove_self_loops_and_parallel_edges) {
    size_t n = G_in.n;
    sequence<tuple<NodeId, NodeId, EdgeTy>> edgelist(G_in.m);
    parallel_for(0, n, [&](NodeId u) {
      parallel_for(G_in.offset[u], G_in.offset[u + 1], [&](EdgeId i) {
        edgelist[i] = make_tuple(u, G_in.E[i].v, G_in.E[i].w);
      });
    });
    sort_inplace(make_slice(edgelist));
    // remove self loops and parallel edges
    auto pred = delayed_seq<bool>(G_in.m, [&](EdgeId i) {
      if (get<0>(edgelist[i]) == get<1>(edgelist[i])) {
        // self loop
        return false;
      } else if (i != 0 && get<0>(edgelist[i]) == get<0>(edgelist[i - 1]) &&
                 get<1>(edgelist[i]) == get<1>(edgelist[i - 1])) {
        // parallel edge
        return false;
      } else {
        return true;
      }
    });
    edgelist = pack(make_slice(edgelist), pred);

    // generate the graph from edgelist
    size_t m = edgelist.size();
    G.n = n;
    G.m = m;
    G.offset = sequence<EdgeId>(n + 1, m);
    G.E = sequence<Edge>(m);
    parallel_for(0, m, [&](EdgeId i) {
      G.E[i] = Edge(get<1>(edgelist[i]), get<2>(edgelist[i]));
      if (i == 0 || get<0>(edgelist[i]) != get<0>(edgelist[i - 1])) {
        G.offset[get<0>(edgelist[i])] = i;
      }
    });
    scan_inclusive_inplace(make_slice(G.offset.rbegin(), G.offset.rend()),
                           minm<EdgeId>());
  } else {
    G.n = G_in.n;
    G.m = G_in.m;
    G.offset = G_in.offset;
    G.E = G_in.E;
  }
  createInEdgesFromOutEdges();
  return;
}

void PCH::pruneNeighbors(NodeId s) {
  if (degree_ordering) {
    if (info[s].in_degree + info[s].out_degree >= g_degree_bound) {
      unsigned long index = newly_inserted_out_edges_hash.first_index(s);
      while (true) {
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_out_edges_hash.H[index].valid) {
          if (newly_inserted_out_edges_hash.H[index].key == s &&
              !vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                       .value.first]) {
            NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
            EdgeTy d = newly_inserted_out_edges_hash.H[index].value.second;
            if (!newly_inserted_out_edges_hash.H[index].settled) {
              newly_inserted_out_edges_hash.H[index].settled = true;
              pair<NodeId, NodeId> nw = make_pair(s, v);
              EdgeTy tentative_dist = vertices_around_hash.find(nw);
              if (tentative_dist < d) {
                newly_inserted_out_edges_hash.H[index].valid = false;
                write_add(&info[s].out_degree, -1);
              }
            }
          }
        }
        index = newly_inserted_out_edges_hash.next_index(index);
      }
      return;
    }
  }
  unordered_map<NodeId, EdgeTy> dist;           // NodeId, distance
  unordered_set<NodeId> neighborsMap;           // NodeId, valid
  unordered_set<NodeId> newSurrondingVertices;  // NodeId, valid
  using T = pair<EdgeTy, NodeId>;  // distance, NodeId, hop_distance
  priority_queue<T, vector<T>, greater<T>> pq;
  dist[s] = 0;
  EdgeTy bound = 0;

  for (size_t i = G.offset[s]; i < G.offset[s + 1]; i++) {
    NodeId v = G.E[i].v;
    if (vertices_contracted[v]) continue;
    dist[v] = G.E[i].w;
    pq.push({dist[v], v});
    if (!early_stop) {
      neighborsMap.insert(v);
      newSurrondingVertices.insert(v);
      if (bound < dist[v]) bound = dist[v];
    }
  }
  unsigned long index = newly_inserted_out_edges_hash.first_index(s);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
    if (newly_inserted_out_edges_hash.H[index].valid) {
      if (newly_inserted_out_edges_hash.H[index].key == s &&
          !vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                   .value.first]) {
        NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
        dist[v] = newly_inserted_out_edges_hash.H[index].value.second;
        pq.push({dist[v], v});
        if (!early_stop) {
          neighborsMap.insert(v);
          newSurrondingVertices.insert(v);
        } else if (!newly_inserted_out_edges_hash.H[index].settled) {
          // newly_inserted_out_edges_hash.H[index].settled = true;
          pair<NodeId, NodeId> nw = make_pair(s, v);
          EdgeTy tentative_dist = vertices_around_hash.find(nw);
          if (tentative_dist < dist[v]) {
            newly_inserted_out_edges_hash.H[index].settled = true;
            newly_inserted_out_edges_hash.H[index].valid = false;
            write_add(&info[s].out_degree, -1);
            // TODO: sync to in_degree
          } else {
            if (bound < dist[v]) bound = dist[v];
            neighborsMap.insert(v);
            newSurrondingVertices.insert(v);
          }
        }
      }
    }
    index = newly_inserted_out_edges_hash.next_index(index);
  }

  NodeId itr = 0;
  while ((!neighborsMap.empty() || !newSurrondingVertices.empty()) &&
         itr < g_max_pop_count && !pq.empty()) {
    auto [d, u] = pq.top();
    if (d > bound) break;
    pq.pop();
    if (dist[u] < d) continue;
    itr++;
    if (degree_ordering) {
      if (info[u].in_degree + info[u].out_degree >= g_degree_bound) {
        continue;
      }
    }
    if (newSurrondingVertices.find(u) != newSurrondingVertices.end()) {
      newSurrondingVertices.erase(u);
      pair<NodeId, NodeId> nw = make_pair(s, u);
      if (vertices_around_hash.find(nw) > d) vertices_around_hash.insert(nw, d);
    }
    bool newly_inserted = (neighborsMap.find(u) != neighborsMap.end());
    if (newly_inserted) {
      neighborsMap.erase(u);
    }
    for (NodeId i = G.offset[u]; i < G.offset[u + 1]; ++i) {
      NodeId v = G.E[i].v;
      if (s == v || vertices_contracted[v]) continue;
      EdgeTy ne = d + G.E[i].w;
      if (dist.find(v) == dist.end() || dist[v] > ne) {
        dist[v] = ne;
        pq.push({ne, v});
        if (newly_inserted) {
          newSurrondingVertices.insert(v);
          if (bound < ne) bound = ne;
        }
      }
    }
    unsigned long index = newly_inserted_out_edges_hash.first_index(s);
    while (true) {
      if (newly_inserted_out_edges_hash.H[index].valid) {
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_out_edges_hash.H[index].key == u &&
            !vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                     .value.first]) {
          NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
          if (s == v) {
            index = newly_inserted_out_edges_hash.next_index(index);
            continue;
          }
          EdgeTy ne = d + newly_inserted_out_edges_hash.H[index].value.second;
          if (dist.find(v) == dist.end() || dist[v] > ne) {
            dist[v] = ne;
            pq.push({ne, v});
            if (newly_inserted) {
              newSurrondingVertices.insert(v);
              if (bound < ne) bound = ne;
            }
          }
        }
      }
      index = newly_inserted_out_edges_hash.next_index(index);
    }
  }
  index = newly_inserted_out_edges_hash.first_index(s);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
    if (newly_inserted_out_edges_hash.H[index].valid) {
      if (newly_inserted_out_edges_hash.H[index].key == s &&
          !vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                   .value.first]) {
        NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
        if (!newly_inserted_out_edges_hash.H[index].settled) {
          newly_inserted_out_edges_hash.H[index].settled = true;
          assert(newly_inserted_out_edges_hash.H[index].value.second >=
                 dist[v]);
          if (newly_inserted_out_edges_hash.H[index].value.second != dist[v]) {
            newly_inserted_out_edges_hash.H[index].valid = false;
            write_add(&info[s].out_degree, -1);
          }
        }
      }
    }
    index = newly_inserted_out_edges_hash.next_index(index);
  }
}

bool PCH::calScore() {
  parallel_for(0, vertices_need_score_update.size(), [&](size_t i) {
    NodeId u = vertices_need_score_update[i];
    if (degree_ordering) {
      if (info[u].in_degree + info[u].out_degree >= g_degree_bound) {
        return;
      }
    }
    vector<NodeId> in_idx(info[u].in_degree);
    vector<NodeId> in_hop(info[u].in_degree);
    vector<EdgeTy> in_wgh(info[u].in_degree);
    vector<NodeId> out_idx(info[u].out_degree);
    vector<NodeId> out_hop(info[u].out_degree);
    vector<EdgeTy> out_wgh(info[u].out_degree);
    transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
    transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
    for (NodeId k1 = 0; k1 < info[u].in_degree; ++k1) {
      for (NodeId k2 = 0; k2 < info[u].out_degree; ++k2) {
        if (in_idx[k1] != out_idx[k2]) {
          pair<NodeId, NodeId> nw = make_pair(in_idx[k1], out_idx[k2]);
          if (vertices_around_hash.find(nw) > in_wgh[k1] + out_wgh[k2]) {
            vertices_around_hash.insert(nw, in_wgh[k1] + out_wgh[k2]);
          }
        }
      }
    }
  });
  parallel_for(0, vertices_need_score_update.size(), [&](size_t i) {
    NodeId u = vertices_need_score_update[i];
    if (degree_ordering) {
      if (info[u].in_degree + info[u].out_degree >= g_degree_bound) {
        info[u].edge_diff = std::numeric_limits<float>::max();
        return;
      }
    }
    if (info[u].in_degree == 0 || info[u].out_degree == 0) {
      info[u].edge_diff = 0;
    } else {
      vector<NodeId> in_idx(info[u].in_degree);
      vector<NodeId> in_hop(info[u].in_degree);
      vector<EdgeTy> in_wgh(info[u].in_degree);
      vector<NodeId> out_idx(info[u].out_degree);
      vector<NodeId> out_hop(info[u].out_degree);
      vector<EdgeTy> out_wgh(info[u].out_degree);
      transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
      transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
      EdgeId removed_arc_num = 1 + info[u].in_degree + info[u].out_degree,
             added_arc_num = 0;
      NodeId removed_hop_count = 1, added_hop_count = 0;
      for (NodeId k1 = 0; k1 < info[u].in_degree; ++k1)
        removed_hop_count += in_hop[k1];
      for (NodeId k1 = 0; k1 < info[u].out_degree; ++k1)
        removed_hop_count += out_hop[k1];
      for (NodeId k1 = 0; k1 < info[u].in_degree; ++k1) {
        for (NodeId k2 = 0; k2 < info[u].out_degree; ++k2) {
          if (in_idx[k1] != out_idx[k2]) {
            pair<NodeId, NodeId> nw = make_pair(in_idx[k1], out_idx[k2]);
            assert(vertices_around_hash.find(nw) <= in_wgh[k1] + out_wgh[k2]);
            if (vertices_around_hash.find(nw) == in_wgh[k1] + out_wgh[k2]) {
              ++added_arc_num;
              added_hop_count += in_hop[k1];
              added_hop_count += out_hop[k2];
            }
          }
        }
      }
      // info[u].edge_diff = 1 + 1000 * G.level[u] +
      //                (1000 * added_arc_num) / removed_arc_num +
      //                (1000 * added_hop_count) / removed_hop_count;
      info[u].edge_diff =
          1 + 500 * G.level[u] +
          1000 * (double)(added_arc_num + 0.2 * info[u].removed_neighbor_num) /
              removed_arc_num +
          1000 * (double)added_hop_count / removed_hop_count +
          10 * (removed_arc_num);
    }
  });
  return true;
}

void PCH::calDisToVerticesAround() {
  // NodeId n = vertices_need_score_update.size();
  // parallel_for(0, n, [&](NodeId i) {
  // pruneNeighbors(vertices_need_score_update[i]); });
  NodeId n = overlay_vertices.size();
  parallel_for(0, n, [&](NodeId i) { pruneNeighbors(overlay_vertices[i]); });
  size_t n2 = bag.pack_into(make_slice(neighbors));
  assert(n2 < neighbors.size());
  parallel_for(0, n2, [&](size_t i) {
    NodeId idx_left = neighbors[i] >> 32;
    NodeId idx_right = neighbors[i] & ((1ull << 32) - 1);
    assert(idx_left != UINT_MAX && idx_right != UINT_MAX);
    // NodeId idx_right = neighbors[i] & (numeric_limits<size_t>::max());
    NodeId key_left = newly_inserted_out_edges_hash.H[idx_left].key;
    NodeId key_right = newly_inserted_in_edges_hash.H[idx_right].key;
    assert(key_left == newly_inserted_in_edges_hash.H[idx_right].value.first);
    assert(key_right == newly_inserted_out_edges_hash.H[idx_left].value.first);
    assert(newly_inserted_out_edges_hash.H[idx_left].value.second ==
           newly_inserted_in_edges_hash.H[idx_right].value.second);
    if (!newly_inserted_out_edges_hash.H[idx_left].valid) {
      if (CAS(&newly_inserted_in_edges_hash.H[idx_right].valid, true, false)) {
        write_add(&info[key_right].in_degree, -1);
      }
    }
  });
  return;
}

template <typename F>
bool PCH::iterateHashTable(NodeId u, F f) {
  unsigned long index = newly_inserted_out_edges_hash.first_index(u);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].valid) {
      if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
      if (newly_inserted_out_edges_hash.H[index].key == u) {
        NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
        if (!f(u, v)) return false;
      }
    }
    index = newly_inserted_out_edges_hash.next_index(index);
  }
  index = newly_inserted_in_edges_hash.first_index(u);
  while (true) {
    if (newly_inserted_in_edges_hash.H[index].valid) {
      if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX) break;
      if (newly_inserted_in_edges_hash.H[index].key == u) {
        NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
        // cerr << "In: ";
        if (!f(u, v)) return false;
      }
    }
    index = newly_inserted_in_edges_hash.next_index(index);
  }
  return true;
}

bool PCH::eligibleForRemoval(NodeId s) {
  auto canRemove = [&](NodeId u, NodeId v) {
    // cerr << "Judge: " << u << ": " << v;
    if (vertices_contracted[v] ||
        info[v].edge_diff > g_upper_score_bound) {
      // cerr << " contracted" << endl;
      return true;
    }

    if ((info[v].edge_diff < info[u].edge_diff) ||
        (info[v].edge_diff == info[u].edge_diff &&
         ((info[v].in_degree + info[v].out_degree) <
              (info[u].in_degree + info[u].out_degree) ||
          ((info[v].in_degree + info[v].out_degree) ==
               (info[u].in_degree + info[u].out_degree) &&
           priority[v] < priority[u])))) {
      // cerr << " False" << endl;
      return false;
    } else {
      // cerr << " Pass" << endl;
      return true;
    }
  };
  assert(info[s].edge_diff <= g_upper_score_bound);
  // printf("In CSR\n");
  for (size_t i = G.offset[s]; i < G.offset[s + 1]; ++i) {
    // cerr << "Out: ";
    if (!canRemove(s, G.E[i].v)) return false;
  }
  for (size_t i = G.in_offset[s]; i < G.in_offset[s + 1]; ++i) {
    // cerr << "In: ";
    if (!canRemove(s, G.in_E[i].v)) return false;
  }
  // printf("In hash table\n");
  return iterateHashTable(s, canRemove);
}

double PCH::sampleUpperBound() {
  size_t seed = hash32(g_tot_level + 1);
  if (overlay_vertices.size() < EDGEDIFF_SAMPLES)
    return info[overlay_vertices[hash32(seed) % overlay_vertices.size()]]
        .edge_diff;
  double sample_edge_diff[EDGEDIFF_SAMPLES + 1];
  for (size_t i = 0; i <= EDGEDIFF_SAMPLES; i++) {
    NodeId v = overlay_vertices[hash32(seed + i) % overlay_vertices.size()];
    assert(!vertices_contracted[v]);
    sample_edge_diff[i] = info[v].edge_diff;
  }
  sort(sample_edge_diff, sample_edge_diff + EDGEDIFF_SAMPLES + 1);
  if (g_sample_bound == 1.0) return std::numeric_limits<double>::max();
  int id = EDGEDIFF_SAMPLES * g_sample_bound;
  return sample_edge_diff[id];
}

bool PCH::insertHelper(NodeId left, NodeId right, EdgeTy len, NodeId hop) {
  assert(left != right);
  assert(hop > 1);
  pair<NodeId, NodeId> nw = make_pair(left, right);
  EdgeTy tentative_dist = vertices_around_hash.find(nw);
  if (tentative_dist < len) return false;
  vertices_settled[left] = false;
  vertices_settled[right] = false;
  bool replaceFlag1 = false;
  bool replaceFlag2 = false;

  for (size_t k = G.offset[left]; k < G.offset[left + 1]; ++k) {
    if (G.E[k].v == right) {
      if (write_min(&G.E[k].w, len, std::less<EdgeTy>())) {
        G.E[k].hop = hop;
      }
      replaceFlag1 = true;
    }
  }
  for (size_t k = G.in_offset[right]; k < G.in_offset[right + 1]; ++k) {
    if (G.in_E[k].v == left) {
      if (write_min(&G.in_E[k].w, len, std::less<EdgeTy>())) {
        G.in_E[k].hop = hop;
      }
      replaceFlag2 = true;
    }
  }
  assert(replaceFlag1 == replaceFlag2);
  if (replaceFlag1) return replaceFlag1;
  pair<bool, NodeId> idx_left =
      newly_inserted_out_edges_hash.insert(left, make_pair(right, len), hop);
  if (idx_left.first) {
    write_add(&info[left].out_degree, 1);
  }
  pair<bool, NodeId> idx_right =
      newly_inserted_in_edges_hash.insert(right, make_pair(left, len), hop);
  if (idx_right.first) {
    write_add(&info[right].in_degree, 1);
  }

  assert(idx_left.second != UINT_MAX && idx_right.second != UINT_MAX);
  bag.insert(((uint64_t)idx_left.second << 32) | idx_right.second);
  return idx_right.first;
}

void PCH::transfer_neighbors(NodeId u, vector<NodeId> &idx, vector<NodeId> &hop, vector<EdgeTy> &wgh,
                             bool forward) {
  const sequence<EdgeId> &offset = forward ? G.offset : G.in_offset;
  const sequence<Edge> &E = forward ? G.E : G.in_E;
  const hash_map<NodeId, NodeId, EdgeTy> &edges_hash_map =
      forward ? newly_inserted_out_edges_hash : newly_inserted_in_edges_hash;
  NodeId degree = forward ? info[u].out_degree : info[u].in_degree;

  NodeId edge_num = 0;
  for (size_t j = offset[u]; j < offset[u + 1]; ++j) {
    if (!vertices_contracted[E[j].v]) {
      idx[edge_num] = E[j].v;
      wgh[edge_num] = E[j].w;
      hop[edge_num] = E[j].hop;
      edge_num++;
    }
  }
  if (edge_num < degree) {
    unsigned long index = edges_hash_map.first_index(u);
    while (true) {
      if (edges_hash_map.H[index].valid) {
        if (edges_hash_map.H[index].key == UINT_N_MAX) break;
        if (edges_hash_map.H[index].key == u &&
            !vertices_contracted[edges_hash_map.H[index].value.first]) {
          idx[edge_num] = edges_hash_map.H[index].value.first;
          wgh[edge_num] = edges_hash_map.H[index].value.second;
          hop[edge_num] = edges_hash_map.H[index].hop;
          edge_num++;
        }
      }
      index = edges_hash_map.next_index(index);
    }
  }
  if (edge_num != degree) {
    printf("FATAL ERROR:\n");
    printf("u: %u, degree: %u, edge_num: %u %s\n", u, degree, edge_num,
           forward ? "outgoing" : "incoming");
    printf("Neighbors of out CSR: ");
    for (size_t j = G.offset[u]; j < G.offset[u + 1]; ++j) {
      printf("(%u,%d,%d) ", G.E[j].v, vertices_contracted[G.E[j].v],
             vertices_contracted[G.E[j].v]);
    }
    puts("");
    printf("Neighbors of in CSR: ");
    for (size_t j = G.in_offset[u]; j < G.in_offset[u + 1]; ++j) {
      printf("(%u,%d,%d) ", G.in_E[j].v, vertices_contracted[G.in_E[j].v],
             vertices_contracted[G.in_E[j].v]);
    }
    puts("");
    printf("Neighbors of out edge hash table: ");
    unsigned long index = newly_inserted_out_edges_hash.first_index(u);
    while (true) {
      if (newly_inserted_out_edges_hash.H[index].valid) {
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_out_edges_hash.H[index].key == u) {
          printf("(%u,%d,%d) ",
                 newly_inserted_out_edges_hash.H[index].value.first,
                 vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                         .value.first],
                 vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                         .value.first]);
        }
      }
      index = newly_inserted_out_edges_hash.next_index(index);
    }
    puts("");
    printf("Neighbors of in edge hash table: ");
    index = newly_inserted_in_edges_hash.first_index(u);
    while (true) {
      if (newly_inserted_in_edges_hash.H[index].valid) {
        if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_in_edges_hash.H[index].key == u) {
          printf("(%u,%d,%d) ",
                 newly_inserted_in_edges_hash.H[index].value.first,
                 vertices_contracted[newly_inserted_in_edges_hash.H[index]
                                         .value.first],
                 vertices_contracted[newly_inserted_in_edges_hash.H[index]
                                         .value.first]);
        }
      }
      index = newly_inserted_in_edges_hash.next_index(index);
    }
    puts("");
    NodeId v = E[offset[u] + 0].v;
    cerr << forward << " " << vertices_contracted[u] << " "
         << vertices_contracted[v] << " " << vertices_contracted[u] << " "
         << vertices_contracted[v] << endl;
    cerr << u << " " << v << " " << edge_num << " " << degree << endl;
    cerr << info[u].edge_diff << " " << info[v].edge_diff << " "
         << info[u].out_degree << " " << info[u].in_degree << " "
         << info[v].out_degree + info[v].in_degree << " " << priority[u] << " "
         << priority[v] << endl;
    assert(edge_num == degree);
  }
  return;
}

void PCH::transferHashtoCSR() {
  GC.offset[GC.n] = 0;
  GC.in_offset[GC.n] = 0;
  parallel_for(0, G.n, [&](size_t i) {
    if (!vertices_contracted[i]) {
      GC.offset[i] = info[i].out_degree;
      GC.in_offset[i] = info[i].in_degree;
    } else {
      GC.offset[i] = 0;
      GC.in_offset[i] = 0;
    }
  });

  GC.m = scan_inplace(GC.offset);
  GC.rm = scan_inplace(GC.in_offset);
  GC.offset[GC.n] = GC.m;
  GC.in_offset[GC.n] = GC.rm;
  GC.E = sequence<Edge>::uninitialized(GC.m);
  GC.in_E = sequence<Edge>::uninitialized(GC.rm);
  parallel_for(0, G.n, [&](size_t i) {
    if (!vertices_contracted[i]) {
      EdgeId bg = GC.offset[i];
      for (EdgeId j = G.offset[i]; j < G.offset[i + 1]; ++j) {
        if (!vertices_contracted[G.E[j].v]) {
          GC.E[bg].v = G.E[j].v;
          GC.E[bg].w = G.E[j].w;
          GC.E[bg].hop = G.E[j].hop;
          bg++;
        }
      }
      bg = clip_from_hash(GC.E, newly_inserted_out_edges_hash, i, bg);
      assert(bg == GC.offset[i + 1]);

      bg = GC.in_offset[i];
      for (EdgeId j = G.in_offset[i]; j < G.in_offset[i + 1]; ++j) {
        if (!vertices_contracted[G.in_E[j].v]) {
          GC.in_E[bg].v = G.in_E[j].v;
          GC.in_E[bg].w = G.in_E[j].w;
          GC.in_E[bg].hop = G.in_E[j].hop;
          bg++;
        }
      }
      bg = clip_from_hash(GC.in_E, newly_inserted_in_edges_hash, i, bg);
      assert(bg == GC.in_offset[i + 1]);
    }
  });
  newly_inserted_out_edges_hash.clear_all();
  newly_inserted_in_edges_hash.clear_all();
  swap(GC.m, G.m);
  swap(GC.E, G.E);
  swap(GC.offset, G.offset);
  swap(GC.rm, G.rm);
  swap(GC.in_E, G.in_E);
  swap(GC.in_offset, G.in_offset);
  cout << "Overlay edges number: " << G.m << endl;
  return;
}

void PCH::checkDegree() {
#ifdef DEBUG
  parallel_for(0, overlay_vertices.size(), [&](size_t nd) {
    NodeId u = overlay_vertices[nd];
    if (vertices_contracted[u]) {
      return;
    }
    size_t in_d = 0, out_d = 0;

    // out_degree
    for (size_t j = G.offset[u]; j < G.offset[u + 1]; j++) {
      NodeId v = G.E[j].v;
      if (!vertices_contracted[v]) {
        out_d++;
      }
    }
    unsigned long index = newly_inserted_out_edges_hash.first_index(u);
    while (true) {
      if (newly_inserted_out_edges_hash.H[index].valid) {
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_out_edges_hash.H[index].key == u) {
          NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
          if (!vertices_contracted[v]) {
            out_d++;
          }
        }
      }
      index = newly_inserted_out_edges_hash.next_index(index);
    }

    // in_degree
    for (size_t j = G.in_offset[u]; j < G.in_offset[u + 1]; j++) {
      NodeId v = G.in_E[j].v;
      if (!vertices_contracted[v]) {
        in_d++;
      }
    }
    index = newly_inserted_in_edges_hash.first_index(u);
    while (true) {
      if (newly_inserted_in_edges_hash.H[index].valid) {
        if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX) break;
        if (newly_inserted_in_edges_hash.H[index].key == u) {
          NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
          if (!vertices_contracted[v]) {
            in_d++;
          }
        }
      }
      index = newly_inserted_in_edges_hash.next_index(index);
    }
    if (info[u].in_degree != in_d) {
      printf("vertex %u has incorrect in_degree, should be %lu but is %u\n", u,
             in_d, info[u].in_degree);
      printf("In neighbors of %u: \n", u);
      // in_degree
      printf("in CSR\n");
      for (size_t j = G.in_offset[u]; j < G.in_offset[u + 1]; j++) {
        NodeId v = G.in_E[j].v;
        printf("(%u,%d) ", v, vertices_contracted[v]);
      }
      printf("\nin hash table\n");
      unsigned long index = newly_inserted_in_edges_hash.first_index(u);
      while (true) {
        if (newly_inserted_in_edges_hash.H[index].valid) {
          if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX) break;
          if (newly_inserted_in_edges_hash.H[index].key == u) {
            NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
            printf("(%u,%d) ", v, vertices_contracted[v]);
          }
        }
        index = newly_inserted_in_edges_hash.next_index(index);
      }
      puts("");
    }
    assert(info[u].in_degree == in_d);
    if (info[u].out_degree != out_d) {
      printf("vertex %u has incorrect out_degree, should be %lu but is %u\n",
             u, out_d, info[u].out_degree);
      printf("Out neighbors of %u: \n", u);
      // out_degree
      printf("in CSR\n");
      for (size_t j = G.offset[u]; j < G.offset[u + 1]; j++) {
        NodeId v = G.E[j].v;
        printf("(%u,%d) ", v, vertices_contracted[v]);
      }
      printf("\nin hash table\n");
      unsigned long index = newly_inserted_out_edges_hash.first_index(u);
      while (true) {
        if (newly_inserted_out_edges_hash.H[index].valid) {
          if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX) break;
          if (newly_inserted_out_edges_hash.H[index].key == u) {
            NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
            printf("(%u,%d) ", v, vertices_contracted[v]);
          }
        } else {
          printf("nv(%u,%d) ",
                 newly_inserted_out_edges_hash.H[index].value.first,
                 vertices_contracted[newly_inserted_out_edges_hash.H[index]
                                         .value.first]);
        }
        index = newly_inserted_out_edges_hash.next_index(index);
      }
      puts("");
    }
    assert(info[u].out_degree == out_d);
  });
#endif
}

bool check_full(hash_map2<NodeId, NodeId, EdgeTy> &vertices_around_hash,
                NodeId g_tot_level) {
  int empty = 0;
  for (int v = 0; v < 20; ++v) {
    if (vertices_around_hash
            .H[hash32(hash32(v) + hash32(g_tot_level)) % vertices_around_hash.m]
            .key.first == numeric_limits<NodeId>::max()) {
      empty++;
    }
  }
  return (empty > 5) ? false : true;
}

void PCH::buildContractionHierarchy() {
  parlay::internal::timer t_contract, t_prune, t_clip, t_score, t_reset,
      t_select;
  t_contract.reset(), t_prune.reset(), t_clip.reset(), t_score.reset(),
      t_reset.reset(), t_select.reset();
  vertices_need_score_update = overlay_vertices;
  NodeId tot = 0;
  bool end_status = false;
  vertices_settled = sequence<bool>(G.n, true);
  vertices_contracted = sequence<bool>(G.n, false);
  ofstream ofs("pch.tsv", ios::app);
  while (overlay_vertices.size() != 0) {
    t_prune.start();
    calDisToVerticesAround();
    t_prune.stop();
    if (end_status) break;
    t_clip.start();
    if (tot > 0.01 * newly_inserted_out_edges_hash.m) {
      checkDegree();
      transferHashtoCSR();

      checkDegree();
      tot = 0;
    }

    t_clip.stop();

    t_score.start();
    calScore();
    t_score.stop();
    t_reset.start();

    g_tot_level++;
    g_upper_score_bound = sampleUpperBound();
    if (g_upper_score_bound == INT_MAX) g_upper_score_bound--;
    parallel_for(0, overlay_vertices.size(), [&](NodeId i) {
      vertices_settled[overlay_vertices[i]] = true;
    });
    t_reset.stop();

    t_select.start();
    auto contracting_vertices = filter(overlay_vertices, [&](NodeId u) {
      return info[u].edge_diff <= g_upper_score_bound && eligibleForRemoval(u);
    });
    t_select.stop();
#ifdef DEBUG
    // printf("checking MIS\n");
    size_t num_can_be_contracted = count_if(overlay_vertices, [&](NodeId u) {
      return info[u].edge_diff <= g_upper_score_bound;
    });
    set<int> st;
    float mx_score = 0, mn_score = INT_MAX;
    NodeId mx_in_degree = 0, mn_in_degree = INT_MAX;
    NodeId mx_out_degree = 0, mn_out_degree = INT_MAX;

    // for (size_t i = 0; i < contracting_vertices.size(); i++) {
    //   mn_score = min(mn_score, info[contracting_vertices[i]].edge_diff);
    // }
    // g_upper_score_bound = mn_score;
    // contracting_vertices = filter(overlay_vertices, [&](NodeId u) {
    //   return info[u].edge_diff <= g_upper_score_bound &&
    //   eligibleForRemoval(u);
    // });
    // num_can_be_contracted = count_if(overlay_vertices, [&](NodeId u) {
    //   return info[u].edge_diff <= g_upper_score_bound;
    // });

    printf("num_can_be_contracted: %zu\n", num_can_be_contracted);
    printf("num_contracting: %zu\n", contracting_vertices.size());
    for (size_t i = 0; i < contracting_vertices.size(); i++) {
      st.insert(contracting_vertices[i]);
      mx_score = max(mx_score, info[contracting_vertices[i]].edge_diff);
      mn_score = min(mn_score, info[contracting_vertices[i]].edge_diff);
      mx_in_degree = max(mx_in_degree, info[contracting_vertices[i]].in_degree);
      mn_in_degree = min(mn_in_degree, info[contracting_vertices[i]].in_degree);
      mx_out_degree =
          max(mx_out_degree, info[contracting_vertices[i]].out_degree);
      mn_out_degree =
          min(mn_out_degree, info[contracting_vertices[i]].out_degree);
    }
    printf("bound: %f, max_score: %f, min_score: %f\n", g_upper_score_bound,
           mx_score, mn_score);
    printf("max_in_degree: %u, min_in_degree: %u\n", mx_in_degree,
           mn_in_degree);
    printf("max_out_degree: %u, min_out_degree: %u\n", mx_out_degree,
           mn_out_degree);
    parallel_for(0, contracting_vertices.size(), [&](size_t i) {
      NodeId u = contracting_vertices[i];
      parallel_for(G.offset[u], G.offset[u + 1],
                   [&](size_t j) { assert(!st.count(G.E[j].v)); });
      parallel_for(G.in_offset[u], G.in_offset[u + 1],
                   [&](size_t j) { assert(!st.count(G.in_E[j].v)); });
      // bool flag =
      //     iterateHashTable(u, [&](NodeId u, NodeId v) { return !st.count(v); });
      // assert(flag == true);
    });
    checkDegree();
#endif
    t_contract.start();
    parallel_for(0, contracting_vertices.size(), [&](size_t i) {
      NodeId u = contracting_vertices[i];
      assert(vertices_contracted[u] == false);
      vector<NodeId> in_idx(info[u].in_degree);
      vector<NodeId> in_hop(info[u].in_degree);
      vector<EdgeTy> in_wgh(info[u].in_degree);
      vector<NodeId> out_idx(info[u].out_degree);
      vector<NodeId> out_hop(info[u].out_degree);
      vector<EdgeTy> out_wgh(info[u].out_degree);
      transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
      transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
      for (NodeId k1 = 0; k1 < info[u].in_degree; ++k1) {
        assert(u != in_idx[k1]);
        for (NodeId k2 = 0; k2 < info[u].out_degree; ++k2) {
          assert(u != out_idx[k2]);
          if (in_idx[k1] != out_idx[k2]) {
            pair<NodeId, NodeId> nw = make_pair(in_idx[k1], out_idx[k2]);
            EdgeTy tentative_dist = vertices_around_hash.find(nw);
            if (tentative_dist >= in_wgh[k1] + out_wgh[k2]) {
              insertHelper(in_idx[k1], out_idx[k2], in_wgh[k1] + out_wgh[k2],
                           in_hop[k1] + out_hop[k2]);
            }
          }
        }
      }
      EdgeId irr = 0;
      for (NodeId k1 = 0; k1 < info[u].in_degree; ++k1) {
        assert(u != in_idx[k1]);
        assert(vertices_contracted[in_idx[k1]] == false);

        pair<NodeId, NodeId> nw = make_pair(in_idx[k1], u);
        EdgeTy tentative_dist = vertices_around_hash.find(nw);
        if (tentative_dist >= in_wgh[k1]) {
          backward_edges_in_ch_hash.insert(u, make_pair(in_idx[k1], in_wgh[k1]),
                                           in_hop[k1]);
          write_max(&G.level[in_idx[k1]], G.level[u] + 1, std::less<int>());
          write_add(&info[in_idx[k1]].removed_neighbor_num, 1);
        } else {
          irr++;
        }
        write_add(&info[in_idx[k1]].out_degree, -1);
      }
      info[u].in_degree -= irr;
      irr = 0;
      for (NodeId k1 = 0; k1 < info[u].out_degree; ++k1) {
        assert(u != out_idx[k1]);
        assert(vertices_contracted[out_idx[k1]] == false);
        pair<NodeId, NodeId> nw = make_pair(u, out_idx[k1]);
        EdgeTy tentative_dist = vertices_around_hash.find(nw);
        if (tentative_dist >= out_wgh[k1]) {
          forward_edges_in_ch_hash.insert(
              u, make_pair(out_idx[k1], out_wgh[k1]), out_hop[k1]);
          write_max(&G.level[out_idx[k1]], G.level[u] + 1, std::less<int>());
          write_add(&info[out_idx[k1]].removed_neighbor_num, 1);
        } else {
          irr++;
        }
        write_add(&info[out_idx[k1]].in_degree, -1);
      }
      info[u].out_degree -= irr;
      vertices_contracted[u] = true;
    });
    t_contract.stop();
    t_reset.start();
    checkDegree();
    NodeId overlay_size = overlay_vertices.size();
    overlay_vertices = parlay::filter(overlay_vertices, [&](NodeId v) {
      return vertices_contracted[v] == 0;
    });
    tot +=
        ((double)G.m / overlay_size) * (overlay_size - overlay_vertices.size());
    if (degree_ordering) {
      if (overlay_size - overlay_vertices.size() == 0) {
        end_status = true;
      }
    }
    if (check_full(vertices_around_hash, g_tot_level)) {
      vertices_need_score_update = overlay_vertices;
      vertices_around_hash.clear_all();
      early_stop = false;
    } else {
      vertices_need_score_update = parlay::filter(
          overlay_vertices,
          [&](NodeId v) { return vertices_settled[v] == false; });
      early_stop = true;
    }
    t_reset.stop();
    if (print_detail) {
      cout << "Round " << g_tot_level << "\n";
      ofs << g_tot_level << "\t";
      cout << "Contract " << contracting_vertices.size() << " vertices\n";
      ofs << contracting_vertices.size() << "\t";
      cout << "Update " << vertices_need_score_update.size() << " vertices\n";
      ofs << vertices_need_score_update.size() << "\t";
      cout << "Overlay vertices number: " << overlay_size << "\n";
      ofs << overlay_size << "\t";
      // EdgeId tot_edge = 0;
      // for (NodeId nid = 0; nid < overlay_vertices.size(); nid++) {
      //   tot_edge += info[overlay_vertices[nid]].in_degree;
      //   tot_edge += info[overlay_vertices[nid]].out_degree;
      // }
      // cout << "Overlay edges number: " << tot_edge << endl;
      // ofs << tot_edge << "\t";
      ofs << t_reset.total_time() << '\t' << t_prune.total_time() << '\t'
          << t_contract.total_time() << '\t' << t_clip.total_time() << '\t'
          << t_score.total_time() << '\t' << t_select.total_time() << '\n';
    }
  }

  ofs << t_reset.total_time() << '\t' << t_prune.total_time() << '\t'
      << t_contract.total_time() << '\t' << t_clip.total_time() << '\t'
      << t_score.total_time() << '\t' << t_select.total_time() << '\t';
  ofs.close();
  cout << "t_reset: " << t_reset.total_time() << '\n';
  cout << "t_prune: " << t_prune.total_time() << '\n';
  cout << "t_contract: " << t_contract.total_time() << '\n';
  cout << "t_clip: " << t_clip.total_time() << '\n';
  cout << "t_score: " << t_score.total_time() << '\n';
}
void PCH::reorderByLayerAndCC(sequence<pair<NodeId, NodeId>> &node_and_ids) {
  if (degree_ordering) {
    parallel_for(0, G.n, [&](size_t i) {
      if (!vertices_contracted[i]) G.level[i] = 0;
    });
  }
  sort_inplace(make_slice(node_and_ids), [&](const pair<NodeId, NodeId> &a,
                                             const pair<NodeId, NodeId> &b) {
    if (a.second != b.second) {
      return a.second < b.second;
    } else if (G.level[a.first] != G.level[b.first]) {
      return G.level[a.first] < G.level[b.first];
    }
    return a.first < b.first;
  });
  G.layerOffset = pack_index<NodeId>(
      delayed_seq<bool>(node_and_ids.size() + 1, [&](size_t i) {
        return i == 0 || i == node_and_ids.size() ||
               node_and_ids[i].second != node_and_ids[i - 1].second ||
               G.level[node_and_ids[i].first] !=
                   G.level[node_and_ids[i - 1].first];
      }));
  G.ccOffset =
      pack_index<NodeId>(delayed_seq<bool>(G.layerOffset.size(), [&](size_t i) {
        return i == 0 || i == G.layerOffset.size() - 1 ||
               node_and_ids[G.layerOffset[i]].second !=
                   node_and_ids[G.layerOffset[i - 1]].second;
      }));
  G.ccRank = sequence<NodeId>(G.n);
  parallel_for(0, G.ccOffset.size() - 1, [&](size_t i) {
    parallel_for(G.layerOffset[G.ccOffset[i]], G.layerOffset[G.ccOffset[i + 1]],
                 [&](size_t j) { G.ccRank[j] = i; });
  });
  return;
}

void reorder(PchGraph &G, bool forward) {
  sequence<EdgeId> &offset = forward ? G.offset : G.in_offset;
  sequence<Edge> &E = forward ? G.E : G.in_E;
  const EdgeId &m = forward ? G.m : G.rm;
  sequence<tuple<NodeId, NodeId, EdgeTy>> edgLst(m, {G.n + 1, G.n + 1, G.n});
  parallel_for(0, G.n, [&](size_t i) {
    parallel_for(offset[i], offset[i + 1], [&](size_t j) {
      edgLst[j] = {G.rank[i], G.rank[E[j].v], E[j].w};
    });
  });
  parallel_for(0, G.n + 1, [&](size_t i) { offset[i] = m; });
  auto smaller2 = [](const tuple<NodeId, NodeId, EdgeTy> &a,
                     const tuple<NodeId, NodeId, EdgeTy> &b) {
    if (get<0>(a) == get<0>(b)) return get<1>(a) < get<1>(b);
    return get<0>(a) < get<0>(b);
  };
  sort_inplace(make_slice(edgLst), smaller2);
  E[0].v = get<1>(edgLst[0]);
  E[0].w = get<2>(edgLst[0]);
  offset[get<0>(edgLst[0])] = 0;
  parallel_for(1, m, [&](size_t i) {
    if (get<0>(edgLst[i]) != get<0>(edgLst[i - 1]))
      offset[get<0>(edgLst[i])] = i;
    E[i].v = get<1>(edgLst[i]);
    E[i].w = get<2>(edgLst[i]);
  });
  parlay::scan_inclusive_inplace(
      parlay::make_slice((offset).rbegin(), (offset).rend()),
      parlay::minm<EdgeId>());
}

void PCH::setOrderedLayer() {
  auto cc_ids = get_cc(G_in);
  sequence<pair<NodeId, NodeId>> node_and_ids = tabulate<pair<NodeId, NodeId>>(
      G_in.n, [&](size_t i) { return make_pair(i, cc_ids[i]); });
  sort_inplace(make_slice(node_and_ids),
               [](const pair<NodeId, NodeId> &a,
                  const pair<NodeId, NodeId> &b) { return a.first < b.first; });
  G.layer = 0;
  G.rank = sequence<NodeId>(G.n);
  G.layer = reduce(make_slice(G.level), maxm<NodeId>());

  reorderByLayerAndCC(node_and_ids);

  parallel_for(0, G.n, [&](size_t i) { G.rank[node_and_ids[i].first] = i; });
  node_and_ids.clear();
  sequence<NodeId> rank_reverse(G.n);
  parallel_for(0, G.n, [&](size_t i) { rank_reverse[G.rank[i]] = G.level[i]; });
  parallel_for(0, G.n, [&](size_t i) { G.level[i] = rank_reverse[i]; });
  rank_reverse.clear();
  reorder(G, true);
  reorder(G, false);
}

void PCH::dumpCH() {
  if (degree_ordering) {
    assert(G.n == GC.n);
    assert(GC.offset.size() == G.n + 1);
    assert(GC.in_offset.size() == G.n + 1);
    assert(G.offset[G.n] == G.m);
    assert(G.in_offset[G.n] == G.m);
    GC.offset[GC.n] = 0;
    GC.in_offset[GC.n] = 0;
    parallel_for(0, G.n, [&](size_t i) {
      GC.offset[i] = info[i].out_degree;
      GC.in_offset[i] = info[i].in_degree;
    });
    GC.m = scan_inplace(GC.offset);
    GC.rm = scan_inplace(GC.in_offset);
    GC.offset[GC.n] = GC.m;
    GC.in_offset[GC.n] = GC.rm;
    GC.E = sequence<Edge>::uninitialized(GC.m);
    GC.in_E = sequence<Edge>::uninitialized(GC.rm);
    parallel_for(0, G.n, [&](size_t i) {
      if (!vertices_contracted[i]) {
        EdgeId bg = GC.offset[i];
        for (EdgeId j = G.offset[i]; j < G.offset[i + 1]; ++j) {
          if (!vertices_contracted[G.E[j].v]) {
            GC.E[bg].v = G.E[j].v;
            GC.E[bg].w = G.E[j].w;
            GC.E[bg].hop = G.E[j].hop;
            bg++;
          }
        }
        bg = clip_from_hash(GC.E, newly_inserted_out_edges_hash, i, bg);
        assert(bg == GC.offset[i + 1]);

        bg = GC.in_offset[i];
        for (EdgeId j = G.in_offset[i]; j < G.in_offset[i + 1]; ++j) {
          if (!vertices_contracted[G.in_E[j].v]) {
            GC.in_E[bg].v = G.in_E[j].v;
            GC.in_E[bg].w = G.in_E[j].w;
            GC.in_E[bg].hop = G.in_E[j].hop;
            bg++;
          }
        }
        bg = clip_from_hash(GC.in_E, newly_inserted_in_edges_hash, i, bg);
        assert(bg == GC.in_offset[i + 1]);
      } else {
        EdgeId bg = GC.offset[i];
        bg = clip_from_hash(GC.E, forward_edges_in_ch_hash, i, bg, false);
        assert(bg == GC.offset[i + 1]);
        bg = GC.in_offset[i];
        bg = clip_from_hash(GC.in_E, backward_edges_in_ch_hash, i, bg, false);
        assert(bg == GC.in_offset[i + 1]);
      }
    });
    swap(GC.m, G.m);
    swap(GC.E, G.E);
    swap(GC.offset, G.offset);
    swap(GC.rm, G.rm);
    swap(GC.in_E, G.in_E);
    swap(GC.in_offset, G.in_offset);
    return;
  }
  G.offset[G.n] = 0;
  G.in_offset[G.n] = 0;
  parallel_for(0, G.n, [&](size_t i) {
    G.offset[i] = info[i].out_degree;
    G.in_offset[i] = info[i].in_degree;
  });

  G.m = scan_inplace(G.offset);
  G.rm = scan_inplace(G.in_offset);
  G.offset[G.n] = G.m;
  G.in_offset[G.n] = G.rm;
  G.E = sequence<Edge>(G.m);
  G.in_E = sequence<Edge>(G.rm);
  parallel_for(0, G.n, [&](size_t i) {
    EdgeId bg = G.offset[i];
    bg = clip_from_hash(G.E, forward_edges_in_ch_hash, i, bg, false);
    assert(bg == G.offset[i + 1]);
    bg = G.in_offset[i];
    bg = clip_from_hash(G.in_E, backward_edges_in_ch_hash, i, bg, false);
    assert(bg == G.in_offset[i + 1]);
  });
}

PchGraph PCH::createContractionHierarchy() {
  internal::timer t;
  t.reset();
  t.start();
  preprocess(true);

#ifdef DEBUG
  parallel_for(0, G.n, [&](size_t i) {
    parallel_for(G.offset[i], G.offset[i + 1], [&](size_t j) {
      assert(i != G.E[j].v);
      if (j > G.offset[i]) {
        assert(G.E[j - 1].v < G.E[j].v);
      }
      NodeId u = i, v = G.E[j].v, w = G.E[j].w;
      auto id = lower_bound(G.in_E.begin() + G.in_offset[v],
                            G.in_E.begin() + G.in_offset[v + 1], Edge(u, 0)) -
                G.in_E.begin();
      assert(G.in_E[id].v == u);
      assert(G.in_E[id].w == w);
    });
  });
#endif
  initialize();
  buildContractionHierarchy();
  dumpCH();
  setOrderedLayer();
  t.stop();
  ofstream ofs("pch.tsv", ios::app);
  ofs << "\t" << G.m << "\t" << G.rm << "\t" << G.m + G.rm << "\t" << G.layer
      << "\t" << t.total_time() << '\n';
  ofs.close();
  printf(
      "G.n: %zu, remaining vertices: %zu, G.m: %zu, G.rm: %zu, G.layer: %zu, "
      "total_time: %f\n",
      G.n, overlay_vertices.size(), G.m, G.rm, G.layer, t.total_time());
  return G;
}
