
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

struct PCH {
private:
  NodeId g_max_pop_count;
  NodeId g_tot_level = 0;
  double g_upper_score_bound = 0;
  sequence<NodeId> overlay_vertices;
  sequence<NodeId> vertices_need_score_update;
  sequence<bool> vertices_settled;
  sequence<bool> vertices_contracted;
  sequence<bool> vertices_contracted_tmp;
  sequence<NodeId> removed_neighbor_num;

  sequence<EdgeId> out_degree; // number of out edges
  sequence<EdgeId> in_degree;  // number of in edges
  sequence<size_t> priority;

  sequence<double> edge_diff;

  hash_map<NodeId, NodeId, EdgeTy> forward_edges_in_ch_hash;
  hash_map<NodeId, NodeId, EdgeTy> backward_edges_in_ch_hash;
  hash_map<NodeId, NodeId, EdgeTy> newly_inserted_out_edges_hash;
  hash_map<NodeId, NodeId, EdgeTy> newly_inserted_in_edges_hash;
  hash_map2<NodeId, NodeId, EdgeTy> vertices_around_hash;
  hashbag<uint64_t> bag;
  sequence<uint64_t> neighbors;

  // helper func
  bool check_edge_valid(size_t idx, NodeId v, bool in_csr);
  void transfer_neighbors(NodeId u, NodeId *idx, NodeId *hop, EdgeTy *wgh,
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
  void pruneNeighbors(NodeId s, bool early_stop);
  bool calScore();
  bool insertHelper(NodeId left, NodeId right, EdgeTy len, NodeId hop,
                    sequence<NodeId> &in_degree_tmp,
                    sequence<NodeId> &out_degree_tmp);
  void reorderByLayerAndCC(sequence<pair<NodeId, NodeId>> &node_and_ids);
  void setOrderedLayer();
  void dumpCH();
  void checkDegree();

  template <typename F> bool iterateHashTable(NodeId u, F f);

public:
  const Graph &G_in;
  PchGraph G;
  PchGraph GC;
  PchGraph createContractionHierarchy();
  PCH(Graph &_G_in, int _max_pop_count = 500)
      : g_max_pop_count(_max_pop_count), removed_neighbor_num(_G_in.n),
        out_degree(_G_in.n), in_degree(_G_in.n), priority(_G_in.n),
        edge_diff(_G_in.n), forward_edges_in_ch_hash(2 * _G_in.m),
        backward_edges_in_ch_hash(2 * _G_in.m),
        newly_inserted_out_edges_hash(_G_in.m),
        newly_inserted_in_edges_hash(_G_in.m),
        vertices_around_hash(10 * _G_in.m), bag(_G_in.n), neighbors(_G_in.n),
        G_in(_G_in) {}
};

bool PCH::check_edge_valid([[maybe_unused]] size_t idx, NodeId v, bool in_csr) {
  if (in_csr) {
    if (vertices_contracted[v])
      return false;
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
      if (hashMap.H[index].key == UINT_N_MAX)
        break;
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
  G.level = sequence<NodeId>(G.n, 0);
  GC.offset = sequence<EdgeId>(GC.n + 1);
  GC.in_offset = sequence<EdgeId>(GC.n + 1);
  GC.E = sequence<Edge>(GC.m);
  GC.in_E = sequence<Edge>(GC.m);
  overlay_vertices = sequence<NodeId>(G.n);
  parallel_for(0, G.n, [&](size_t i) {
    overlay_vertices[i] = i;
    out_degree[i] = G.offset[i + 1] - G.offset[i];
    in_degree[i] = G.in_offset[i + 1] - G.in_offset[i];
    priority[i] = i + 1;
    removed_neighbor_num[i] = 0;
  });
  priority = parlay::random_shuffle(priority);
}

void PCH::createInEdgesFromOutEdges() {
  size_t n = G.n;
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
  // G.in_E[0].v = get<1>(edgelist[0]);
  // G.in_E[0].w = get<2>(edgelist[0]);
  // G.in_offset[get<0>(edgelist[0])] = 0;
  // parallel_for(1, G.m, [&](size_t i) {
  //   if (get<0>(edgelist[i]) != get<0>(edgelist[i - 1]))
  //     G.in_offset[get<0>(edgelist[i])] = i;
  //   G.in_E[i].v = get<1>(edgelist[i]);
  //   G.in_E[i].w = get<2>(edgelist[i]);
  // });
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

void PCH::pruneNeighbors(NodeId s, bool early_stop) {
  unordered_map<NodeId, EdgeTy> dist;          // NodeId, distance
  unordered_set<NodeId> neighborsMap;          // NodeId, valid
  unordered_set<NodeId> newSurrondingVertices; // NodeId, valid
  using T = pair<EdgeTy, NodeId>;              // distance, NodeId, hop_distance
  priority_queue<T, vector<T>, greater<T>> pq;
  dist[s] = 0;
  EdgeTy bound = 0;

  for (size_t i = G.offset[s]; i < G.offset[s + 1]; i++) {
    NodeId v = G.E[i].v;
    if (vertices_contracted[v])
      continue;
    dist[v] = G.E[i].w;
    pq.push({dist[v], v});
    if (!early_stop) {
      neighborsMap.insert(v);
      newSurrondingVertices.insert(v);
      if (bound < dist[v])
        bound = dist[v];
    }
  }
  unsigned long index = newly_inserted_out_edges_hash.first_index(s);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
      break;
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
            write_add(&out_degree[s], -1);
            // TODO: sync to in_degree
          } else {
            if (bound < dist[v])
              bound = dist[v];
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
    if (d > bound)
      break;
    pq.pop();
    if (dist[u] < d)
      continue;
    itr++;
    if (newSurrondingVertices.find(u) != newSurrondingVertices.end()) {
      newSurrondingVertices.erase(u);
      pair<NodeId, NodeId> nw = make_pair(s, u);
      if (vertices_around_hash.find(nw) > d)
        vertices_around_hash.insert(nw, d);
    }
    bool newly_inserted = (neighborsMap.find(u) != neighborsMap.end());
    if (newly_inserted) {
      neighborsMap.erase(u);
    }
    for (NodeId i = G.offset[u]; i < G.offset[u + 1]; ++i) {
      NodeId v = G.E[i].v;
      if (s == v || vertices_contracted[v])
        continue;
      EdgeTy ne = d + G.E[i].w;
      if (dist.find(v) == dist.end() || dist[v] > ne) {
        dist[v] = ne;
        pq.push({ne, v});
        if (newly_inserted) {
          newSurrondingVertices.insert(v);
          if (bound < ne)
            bound = ne;
        }
      }
    }
    unsigned long index = newly_inserted_out_edges_hash.first_index(s);
    while (true) {
      if (newly_inserted_out_edges_hash.H[index].valid) {
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
          break;
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
              if (bound < ne)
                bound = ne;
            }
          }
        }
      }
      index = newly_inserted_out_edges_hash.next_index(index);
    }
  }
  index = newly_inserted_out_edges_hash.first_index(s);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
      break;
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
            write_add(&out_degree[s], -1);
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
    NodeId in_idx[in_degree[u]];
    NodeId in_hop[in_degree[u]];
    EdgeTy in_wgh[in_degree[u]];
    NodeId out_idx[out_degree[u]];
    NodeId out_hop[out_degree[u]];
    EdgeTy out_wgh[out_degree[u]];
    transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
    transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
    for (NodeId k1 = 0; k1 < in_degree[u]; ++k1) {
      for (NodeId k2 = 0; k2 < out_degree[u]; ++k2) {
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
    NodeId in_idx[in_degree[u]];
    NodeId in_hop[in_degree[u]];
    EdgeTy in_wgh[in_degree[u]];
    NodeId out_idx[out_degree[u]];
    NodeId out_hop[out_degree[u]];
    EdgeTy out_wgh[out_degree[u]];
    transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
    transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
    EdgeId removed_arc_num = 1 + in_degree[u] + out_degree[u],
           added_arc_num = 0;
    NodeId removed_hop_count = 1, added_hop_count = 0;
    for (NodeId k1 = 0; k1 < in_degree[u]; ++k1)
      removed_hop_count += in_hop[k1];
    for (NodeId k1 = 0; k1 < out_degree[u]; ++k1)
      removed_hop_count += out_hop[k1];
    for (NodeId k1 = 0; k1 < in_degree[u]; ++k1) {
      for (NodeId k2 = 0; k2 < out_degree[u]; ++k2) {
        if (in_idx[k1] != out_idx[k2]) {
          pair<NodeId, NodeId> nw = make_pair(in_idx[k1], out_idx[k2]);
          assert(vertices_around_hash.find(nw) <= in_wgh[k1] + out_wgh[k2]);
          if (vertices_around_hash.find(nw) == in_wgh[k1] + out_wgh[k2]) {
            // if (vertices_around_hash.find(nw) > in_wgh[k1] + out_wgh[k2]) {
            //   vertices_around_hash.insert(nw, in_wgh[k1] + out_wgh[k2]);
            // }
            ++added_arc_num;
            added_hop_count += in_hop[k1];
            added_hop_count += out_hop[k2];
          }
        }
      }
    }
    // edge_diff[u] = 1 + 1000 * G.level[u] +
    //                (1000 * added_arc_num) / removed_arc_num +
    //                (1000 * added_hop_count) / removed_hop_count;
    edge_diff[u] = 1 + 1000 * G.level[u] +
                   1000 *
                       (double)(added_arc_num + 0.2 * removed_neighbor_num[u]) /
                       removed_arc_num +
                   1000 * (double)added_hop_count / removed_hop_count +
                   10 * (removed_arc_num);
  });
  return true;
}

void PCH::calDisToVerticesAround() {
  NodeId n = vertices_need_score_update.size();
  // cerr << "n = " << n << endl;
  parallel_for(0, n, [&](NodeId i) {
    pruneNeighbors(vertices_need_score_update[i], (g_tot_level != 0));
    // vertices_settled[vertices_need_score_update[i]] = true;
  });
  // TODO: sync pruned outedgree to in_degree
  size_t n2 = bag.pack_into(make_slice(neighbors));
  // cerr << "n2 = " << n2 << endl;
  parallel_for(0, n2, [&](size_t i) {
    NodeId idx_left = neighbors[i] >> 32;
    NodeId idx_right = neighbors[i] & ((1ull << 32) - 1);
    // NodeId idx_right = neighbors[i] & (numeric_limits<size_t>::max());
    NodeId key_left = newly_inserted_out_edges_hash.H[idx_left].key;
    NodeId key_right = newly_inserted_in_edges_hash.H[idx_right].key;
    assert(key_left == newly_inserted_in_edges_hash.H[idx_right].value.first);
    assert(key_right == newly_inserted_out_edges_hash.H[idx_left].value.first);
    assert(newly_inserted_out_edges_hash.H[idx_left].value.second ==
           newly_inserted_in_edges_hash.H[idx_right].value.second);
    // if (key_left == 833739 ||  key_right == 833739) {
    //   cerr << "Try to delete " << key_left << " " << key_right << " "
    //        << newly_inserted_out_edges_hash.H[idx_left].valid << " "
    //        << newly_inserted_in_edges_hash.H[idx_right].valid << endl;
    // }
    if (!newly_inserted_out_edges_hash.H[idx_left].valid) {
      if (CAS(&newly_inserted_in_edges_hash.H[idx_right].valid, true, false)) {
        write_add(&in_degree[key_right], -1);
      }
    }
  });
  return;
}

template <typename F> bool PCH::iterateHashTable(NodeId u, F f) {
  unsigned long index = newly_inserted_out_edges_hash.first_index(u);
  while (true) {
    if (newly_inserted_out_edges_hash.H[index].valid) {
      if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
        break;
      if (newly_inserted_out_edges_hash.H[index].key == u) {
        NodeId v = newly_inserted_out_edges_hash.H[index].value.first;
        // cerr << "Out: ";
        if (!f(u, v))
          return false;
      }
    }
    index = newly_inserted_out_edges_hash.next_index(index);
  }
  index = newly_inserted_in_edges_hash.first_index(u);
  while (true) {
    if (newly_inserted_in_edges_hash.H[index].valid) {
      if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX)
        break;
      if (newly_inserted_in_edges_hash.H[index].key == u) {
        NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
        // cerr << "In: ";
        if (!f(u, v))
          return false;
      }
    }
    index = newly_inserted_in_edges_hash.next_index(index);
  }
  return true;
}

bool PCH::eligibleForRemoval(NodeId s) {
  auto canRemove = [&](NodeId u, NodeId v) {
    // cerr << "Judge: " << u << ": " << v;
    if (vertices_contracted[v] || edge_diff[v] > g_upper_score_bound) {
      // cerr << " contracted" << endl;
      return true;
    }

    // if (edge_diff[v] <= g_upper_score_bound) {
    //   if (edge_diff[v] == edge_diff[u]) {
    //     return priority[v] < priority[u];
    //   }
    //   return edge_diff[v] > edge_diff[u];
    // } else {
    //   return true;
    // }

    if ((edge_diff[v] < edge_diff[u]) ||
        (edge_diff[v] == edge_diff[u] &&
         ((in_degree[v] + out_degree[v]) < (in_degree[u] + out_degree[u]) ||
          ((in_degree[v] + out_degree[v]) == (in_degree[u] + out_degree[u]) &&
           priority[v] < priority[u])))) {
      // cerr << " False" << endl;
      return false;
    } else {
      // cerr << " Pass" << endl;
      return true;
    }
  };
  assert(edge_diff[s] <= g_upper_score_bound);
  // printf("In CSR\n");
  for (size_t i = G.offset[s]; i < G.offset[s + 1]; ++i) {
    // cerr << "Out: ";
    if (!canRemove(s, G.E[i].v))
      return false;
  }
  for (size_t i = G.in_offset[s]; i < G.in_offset[s + 1]; ++i) {
    // cerr << "In: ";
    if (!canRemove(s, G.in_E[i].v))
      return false;
  }
  // printf("In hash table\n");
  return iterateHashTable(s, canRemove);
}

double PCH::sampleUpperBound() {
  size_t seed = hash32(g_tot_level + 1);
  if (overlay_vertices.size() < EDGEDIFF_SAMPLES)
    return edge_diff[overlay_vertices[hash32(seed) % overlay_vertices.size()]];
  double sample_edge_diff[EDGEDIFF_SAMPLES + 1];
  for (size_t i = 0; i <= EDGEDIFF_SAMPLES; i++) {
    NodeId v = overlay_vertices[hash32(seed + i) % overlay_vertices.size()];
    assert(!vertices_contracted[v]);
    sample_edge_diff[i] = edge_diff[v];
  }
  sort(sample_edge_diff, sample_edge_diff + EDGEDIFF_SAMPLES + 1);
  // cout << sample_edge_diff[0] << "\t" << sample_edge_diff[10] << "\t"
  //      << sample_edge_diff[100] << endl;
  int id = EDGEDIFF_SAMPLES * 0.1;
  return sample_edge_diff[id];
}

bool PCH::insertHelper(NodeId left, NodeId right, EdgeTy len, NodeId hop,
                       sequence<NodeId> &in_degree_tmp,
                       sequence<NodeId> &out_degree_tmp) {
  assert(left != right);
  assert(hop > 1);
  pair<NodeId, NodeId> nw = make_pair(left, right);
  EdgeTy tentative_dist = vertices_around_hash.find(nw);
  if (tentative_dist < len)
    return false;
  // cerr << "Insert " << left << " " << right << " " << len << endl;
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
  if (replaceFlag1)
    return replaceFlag1;
  pair<bool, NodeId> idx_left =
      newly_inserted_out_edges_hash.insert(left, make_pair(right, len), hop);
  if (idx_left.first) {
    write_add(&out_degree_tmp[left], 1);
  }
  pair<bool, NodeId> idx_right =
      newly_inserted_in_edges_hash.insert(right, make_pair(left, len), hop);
  if (idx_right.first) {
    write_add(&in_degree_tmp[right], 1);
  }

  bag.insert(((uint64_t)idx_left.second << 32) | idx_right.second);
  return idx_right.first;
}

void PCH::transfer_neighbors(NodeId u, NodeId *idx, NodeId *hop, EdgeTy *wgh,
                             bool forward) {
  const sequence<EdgeId> &offset = forward ? G.offset : G.in_offset;
  const sequence<Edge> &E = forward ? G.E : G.in_E;
  const hash_map<NodeId, NodeId, EdgeTy> &edges_hash_map =
      forward ? newly_inserted_out_edges_hash : newly_inserted_in_edges_hash;
  NodeId degree = forward ? out_degree[u] : in_degree[u];

  NodeId edge_num = 0;
  for (size_t j = offset[u]; j < offset[u + 1]; ++j) {
    if (!vertices_contracted[E[j].v]) {
      idx[edge_num] = E[j].v;
      wgh[edge_num] = E[j].w;
      hop[edge_num] = E[j].hop;
      edge_num++;
      // if (u == 39 || u == 21 || u == 58)
      //   cerr << u << ": " << E[j].v << " " << forward << endl;
    }
  }
  if (edge_num < degree) {
    unsigned long index = edges_hash_map.first_index(u);
    while (true) {
      if (edges_hash_map.H[index].valid) {
        if (edges_hash_map.H[index].key == UINT_N_MAX)
          break;
        if (edges_hash_map.H[index].key == u &&
            !vertices_contracted[edges_hash_map.H[index].value.first]) {
          idx[edge_num] = edges_hash_map.H[index].value.first;
          wgh[edge_num] = edges_hash_map.H[index].value.second;
          hop[edge_num] = edges_hash_map.H[index].hop;
          edge_num++;
          // if (u == 39 || u == 21 || u == 58)
          //   cerr << u << ": " << edges_hash_map.H[index].value.first << " "
          //        << forward << endl;
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
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
          break;
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
        if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX)
          break;
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
    cerr << edge_diff[u] << " " << edge_diff[v] << " " << out_degree[u] << " "
         << in_degree[u] << " " << out_degree[v] + in_degree[v] << " "
         << priority[u] << " " << priority[v] << endl;
    assert(edge_num == degree);
  }
  return;
}

void PCH::transferHashtoCSR() {
  GC.offset[GC.n] = 0;
  GC.in_offset[GC.n] = 0;
  parallel_for(0, G.n, [&](size_t i) {
    if (!vertices_contracted[i]) {
      GC.offset[i] = out_degree[i];
      GC.in_offset[i] = in_degree[i];
    } else {
      GC.offset[i] = 0;
      GC.in_offset[i] = 0;
    }
  });

  GC.m = scan_inplace(GC.offset);
  GC.rm = scan_inplace(GC.in_offset);
  GC.offset[GC.n] = GC.m;
  GC.in_offset[GC.n] = GC.rm;
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
  // newly_inserted_out_edges_hash.m = max((size_t)GC.n / 500, (size_t)(GC.m));
  // newly_inserted_in_edges_hash.m = max((size_t)GC.n / 500, (size_t)(GC.rm));
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
        if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
          break;
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
        if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX)
          break;
        if (newly_inserted_in_edges_hash.H[index].key == u) {
          NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
          if (!vertices_contracted[v]) {
            in_d++;
          }
        }
      }
      index = newly_inserted_in_edges_hash.next_index(index);
    }
    if (in_degree[u] != in_d) {
      printf("vertex %u has incorrect in_degree, should be %lu but is %zu\n", u,
             in_d, in_degree[u]);
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
          if (newly_inserted_in_edges_hash.H[index].key == UINT_N_MAX)
            break;
          if (newly_inserted_in_edges_hash.H[index].key == u) {
            NodeId v = newly_inserted_in_edges_hash.H[index].value.first;
            printf("(%u,%d) ", v, vertices_contracted[v]);
          }
        }
        index = newly_inserted_in_edges_hash.next_index(index);
      }
      puts("");
    }
    assert(in_degree[u] == in_d);
    if (out_degree[u] != out_d) {
      printf("vertex %u has incorrect out_degree, should be %lu but is %zu\n",
             u, out_d, out_degree[u]);
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
          if (newly_inserted_out_edges_hash.H[index].key == UINT_N_MAX)
            break;
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
    assert(out_degree[u] == out_d);
  });
#endif
}

void PCH::buildContractionHierarchy() {
  parlay::internal::timer t_contract, t_prune, t_clip, t_score, t_reset;
  t_contract.reset(), t_prune.reset(), t_clip.reset(), t_score.reset(),
      t_reset.reset();
  // bool continue_contraction_process = true;
  vertices_need_score_update = overlay_vertices;
  sequence<NodeId> in_degree_tmp(G.n);
  sequence<NodeId> out_degree_tmp(G.n);
  vertices_settled = sequence<bool>(G.n, true);
  vertices_contracted = sequence<bool>(G.n, false);
  vertices_contracted_tmp = sequence<bool>(G.n, false);
  NodeId tot = 0;
  // int itr = 0;
  cout << "Overlay vertices number: " << overlay_vertices.size() << endl;
  while (overlay_vertices.size() != 0) {
    // cerr << itr++ << endl;
    t_prune.start();
    calDisToVerticesAround();
    t_prune.stop();
    // cerr << "calDisToVerticesAround" << endl;

    t_clip.start();
    if (tot > 0.01 * newly_inserted_out_edges_hash.m) {
      // printf("before clip\n");
      checkDegree();
      transferHashtoCSR();

      checkDegree();
      // printf("after clip\n");
      tot = 0;
    }

    t_clip.stop();
    // cerr << "clip" << endl;

    t_score.start();
    calScore();
    t_score.stop();
    // cerr << "calScore" << endl;

    t_contract.start();
    g_tot_level++;
    g_upper_score_bound = sampleUpperBound();
    parallel_for(0, overlay_vertices.size(), [&](NodeId i) {
      NodeId u = overlay_vertices[i];
      in_degree_tmp[u] = in_degree[u];
      out_degree_tmp[u] = out_degree[u];
      vertices_contracted_tmp[u] = vertices_contracted[u];
      vertices_settled[u] = true;
    });
    // cerr << "sampleAndDuplicate " << g_upper_score_bound << endl;
    // printf("Start Contracting\n");

#ifdef DEBUG
    // printf("checking MIS\n");
    auto contracting_vertices = filter(overlay_vertices, [&](NodeId u) {
      return edge_diff[u] <= g_upper_score_bound && eligibleForRemoval(u);
    });
    size_t num_can_be_contracted = count_if(overlay_vertices, [&](NodeId u) {
      return edge_diff[u] <= g_upper_score_bound;
    });
    set<int> st;
    double mx_score = 0, mn_score = INT_MAX;
    EdgeId mx_in_degree = 0, mn_in_degree = INT_MAX;
    EdgeId mx_out_degree = 0, mn_out_degree = INT_MAX;

    // for (size_t i = 0; i < contracting_vertices.size(); i++) {
    //   mn_score = min(mn_score, edge_diff[contracting_vertices[i]]);
    // }
    // g_upper_score_bound = mn_score;
    // contracting_vertices = filter(overlay_vertices, [&](NodeId u) {
    //   return edge_diff[u] <= g_upper_score_bound && eligibleForRemoval(u);
    // });
    // num_can_be_contracted = count_if(overlay_vertices, [&](NodeId u) {
    //   return edge_diff[u] <= g_upper_score_bound;
    // });

    printf("num_can_be_contracted: %zu\n", num_can_be_contracted);
    printf("num_contracting: %zu\n", contracting_vertices.size());
    for (size_t i = 0; i < contracting_vertices.size(); i++) {
      st.insert(contracting_vertices[i]);
      mx_score = max(mx_score, edge_diff[contracting_vertices[i]]);
      mn_score = min(mn_score, edge_diff[contracting_vertices[i]]);
      mx_in_degree = max(mx_in_degree, in_degree[contracting_vertices[i]]);
      mn_in_degree = min(mn_in_degree, in_degree[contracting_vertices[i]]);
      mx_out_degree = max(mx_out_degree, out_degree[contracting_vertices[i]]);
      mn_out_degree = min(mn_out_degree, out_degree[contracting_vertices[i]]);
    }
    printf("bound: %f, max_score: %f, min_score: %f\n", g_upper_score_bound,
           mx_score, mn_score);
    printf("max_in_degree: %lu, min_in_degree: %lu\n", mx_in_degree,
           mn_in_degree);
    printf("max_out_degree: %lu, min_out_degree: %lu\n", mx_out_degree,
           mn_out_degree);
    parallel_for(0, contracting_vertices.size(), [&](size_t i) {
      NodeId u = contracting_vertices[i];
      parallel_for(G.offset[u], G.offset[u + 1],
                   [&](size_t j) { assert(!st.count(G.E[j].v)); });
      parallel_for(G.in_offset[u], G.in_offset[u + 1],
                   [&](size_t j) { assert(!st.count(G.in_E[j].v)); });
      bool flag =
          iterateHashTable(u, [&](NodeId u, NodeId v) { return !st.count(v); });
      assert(flag == true);
    });
#endif
    // printf("before contracting\n");
    checkDegree();
    // printf("before contracting checked\n");
    parallel_for(0, overlay_vertices.size(), [&](size_t i) {
      NodeId u = overlay_vertices[i];
      if (edge_diff[u] <= g_upper_score_bound) {
        if (eligibleForRemoval(u)) {
          // printf("Node %u eligible\n", u);
          assert(vertices_contracted_tmp[u] == false);
          NodeId in_idx[in_degree[u]];
          NodeId in_hop[in_degree[u]];
          EdgeTy in_wgh[in_degree[u]];
          NodeId out_idx[out_degree[u]];
          NodeId out_hop[out_degree[u]];
          EdgeTy out_wgh[out_degree[u]];
          transfer_neighbors(u, in_idx, in_hop, in_wgh, false);
          transfer_neighbors(u, out_idx, out_hop, out_wgh, true);
          for (NodeId k1 = 0; k1 < in_degree[u]; ++k1) {
            assert(u != in_idx[k1]);
            for (NodeId k2 = 0; k2 < out_degree[u]; ++k2) {
              assert(u != out_idx[k2]);
              if (in_idx[k1] != out_idx[k2]) {
                pair<NodeId, NodeId> nw = make_pair(in_idx[k1], out_idx[k2]);
                EdgeTy tentative_dist = vertices_around_hash.find(nw);
                if (tentative_dist >= in_wgh[k1] + out_wgh[k2]) {
                  insertHelper(
                      in_idx[k1], out_idx[k2], in_wgh[k1] + out_wgh[k2],
                      in_hop[k1] + out_hop[k2], in_degree_tmp, out_degree_tmp);
                }
              }
            }
          }
          for (NodeId k1 = 0; k1 < in_degree[u]; ++k1) {
            assert(u != in_idx[k1]);
            assert(vertices_contracted[in_idx[k1]] == false);

            pair<NodeId, NodeId> nw = make_pair(in_idx[k1], u);
            EdgeTy tentative_dist = vertices_around_hash.find(nw);
            if (tentative_dist >= in_wgh[k1]) {
              backward_edges_in_ch_hash.insert(
                  u, make_pair(in_idx[k1], in_wgh[k1]), in_hop[k1]);
              write_max(&G.level[in_idx[k1]], G.level[u] + 1, std::less<int>());
              write_add(&removed_neighbor_num[in_idx[k1]], 1);
            } else {
              write_add(&in_degree_tmp[u], -1);
            }
            write_add(&out_degree_tmp[in_idx[k1]], -1);
          }
          for (NodeId k1 = 0; k1 < out_degree[u]; ++k1) {
            assert(u != out_idx[k1]);
            assert(vertices_contracted[out_idx[k1]] == false);
            pair<NodeId, NodeId> nw = make_pair(u, out_idx[k1]);
            EdgeTy tentative_dist = vertices_around_hash.find(nw);
            if (tentative_dist >= out_wgh[k1]) {
              forward_edges_in_ch_hash.insert(
                  u, make_pair(out_idx[k1], out_wgh[k1]), out_hop[k1]);
              write_max(&G.level[out_idx[k1]], G.level[u] + 1, std::less<int>());
              write_add(&removed_neighbor_num[out_idx[k1]], 1);
            } else {
              write_add(&out_degree_tmp[u], -1);
            }
            write_add(&in_degree_tmp[out_idx[k1]], -1);
          }
          vertices_contracted_tmp[u] = true;
        }
      }
    });
    t_contract.stop();
    // cerr << "Contraction" << endl;
    t_reset.start();
    parallel_for(0, overlay_vertices.size(), [&](size_t nd) {
      NodeId i = overlay_vertices[nd];
      in_degree[i] = in_degree_tmp[i];
      out_degree[i] = out_degree_tmp[i];
      vertices_contracted[i] = vertices_contracted_tmp[i];
    });
    // printf("after contracting\n");
    checkDegree();
    // printf("after contracting checked\n");
    sequence<NodeId> overlay_vertices_tmp =
        parlay::filter(overlay_vertices,
                       [&](NodeId v) { return vertices_contracted[v] == 0; });
    tot += ((double)G.m / overlay_vertices.size()) *
           (overlay_vertices.size() - overlay_vertices_tmp.size());
    swap(overlay_vertices_tmp, overlay_vertices);
    vertices_need_score_update =
        parlay::filter(overlay_vertices,
                       [&](NodeId v) { return vertices_settled[v] == false; });
    t_reset.stop();
    cout << "Overlay vertices number: " << overlay_vertices.size() << endl;
    // cerr << "reset" << endl;
  }
  
  cout<<t_reset.total_time()<<"\t";
  cout<<t_prune.total_time()<<"\t";
  cout<<t_contract.total_time()<<"\t";
  cout<<t_clip.total_time()<<"\t";
  cout<<t_score.total_time()<<endl;
}
void PCH::reorderByLayerAndCC(sequence<pair<NodeId, NodeId>> &node_and_ids) {
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
  G.layerOffset[G.layerOffset.size() - 1] = G.n;
  G.ccOffset = pack_index<NodeId>(
      delayed_seq<bool>(G.layerOffset.size() + 1, [&](size_t i) {
        return i == 0 || i == G.layerOffset.size() ||
               node_and_ids[G.layerOffset[i]].second !=
                   node_and_ids[G.layerOffset[i - 1]].second;
      }));
  G.ccOffset[G.ccOffset.size() - 1] = G.layerOffset.size() - 1;
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
    if (get<0>(a) == get<0>(b))
      return get<1>(a) < get<1>(b);
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
  sequence<std::pair<size_t, NodeId>> layerMap(G.n);
  G.rank = sequence<NodeId>(G.n);
  G.layer = reduce(make_slice(G.level), maxm<NodeId>());

  reorderByLayerAndCC(node_and_ids);

  parallel_for(0, G.n, [&](size_t i) { G.rank[node_and_ids[i].first] = i; });
  // TODO(xiaojun): use tabulate
  // G.level = tabulate(G.n, [&](size_t i) {
  //   return G.level[G.order[i]];
  // });
  sequence<NodeId> rank_reverse(G.n);
  parallel_for(0, G.n, [&](size_t i) { rank_reverse[G.rank[i]] = G.level[i]; });
  parallel_for(0, G.n, [&](size_t i) { G.level[i] = rank_reverse[i]; });
  reorder(G, true);
  reorder(G, false);
}

void PCH::dumpCH() {
  G.offset[G.n] = 0;
  G.in_offset[G.n] = 0;
  parallel_for(0, G.n, [&](size_t i) {
    G.offset[i] = out_degree[i];
    G.in_offset[i] = in_degree[i];
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
  cout << G.n << "\t" << G.m << "\t" << G.rm  << "\t" << G.layer << "\t" << t.total_time() << endl;
  return G;
}
