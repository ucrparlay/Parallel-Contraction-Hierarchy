#include "graph.hpp"

struct PchQuery {
  size_t n;
  const PchGraph &GC;
  const Graph &G;
  sequence<EdgeTy> dist;
  sequence<uint32_t> timestamp;
  uint32_t current_timestamp;

  // for testing
  sequence<EdgeId> in_offset_incoming;
  sequence<Edge> in_E_incoming;

  PchQuery(const PchGraph &_GC, const Graph &_G) : GC(_GC), G(_G) {
    assert(G.n == GC.n);
    n = G.n;
    dist = sequence<EdgeTy>(n * 2);
    timestamp = sequence<uint32_t>(n * 2);
    current_timestamp = 0;
  }

  sequence<EdgeTy> ssspVerifier(NodeId s);
  pair<EdgeTy, int> stVerifier(NodeId s, NodeId t);
  pair<EdgeTy, int> stQuery(NodeId s, NodeId t);
  sequence<EdgeTy> ssspQuery(NodeId s, bool remap, bool in_parallel);

  // for testing
  void make_inverse();
};

void PchQuery::make_inverse() {
  parlay::sequence<std::pair<NodeId, Edge>> edgelist(GC.rm);
  parlay::parallel_for(0, GC.n, [&](NodeId u) {
    parlay::parallel_for(GC.in_offset[u], GC.in_offset[u + 1], [&](EdgeId i) {
      edgelist[i] = std::make_pair(GC.in_E[i].v, Edge(u, GC.in_E[i].w));
    });
  });
  parlay::sort_inplace(
      parlay::make_slice(edgelist),
      [](const std::pair<NodeId, Edge> &a, const std::pair<NodeId, Edge> &b) {
        if (a.first != b.first) {
          return a.first < b.first;
        }
        return a.second.v < b.second.v;
      });
  in_offset_incoming = parlay::sequence<EdgeId>(n + 1, GC.rm);
  in_E_incoming = parlay::sequence<Edge>(GC.rm);
  parlay::parallel_for(0, GC.rm, [&](size_t i) {
    in_E_incoming[i] = edgelist[i].second;
    if (i == 0 || edgelist[i].first != edgelist[i - 1].first) {
      in_offset_incoming[edgelist[i].first] = i;
    }
  });
  parlay::scan_inclusive_inplace(parlay::make_slice(in_offset_incoming.rbegin(),
                                                    in_offset_incoming.rend()),
                                 parlay::minm<EdgeId>());
}

sequence<EdgeTy> PchQuery::ssspVerifier(NodeId s) {
  sequence<EdgeTy> dist(n, DIST_MAX);
  using P = pair<EdgeTy, NodeId>;
  priority_queue<P, vector<P>, greater<P>> pq;
  dist[s] = 0;
  pq.push({dist[s], s});
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();
    if (d != dist[u]) {
      continue;
    }
    for (size_t i = G.offset[u]; i < G.offset[u + 1]; i++) {
      NodeId v = G.E[i].v;
      EdgeTy w = G.E[i].w;
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push({dist[v], v});
      }
    }
  }
  return dist;
}

pair<EdgeTy, int> PchQuery::stVerifier(NodeId s, NodeId t) {
  // TODO: change to bidirectional search
  sequence<EdgeTy> dist(n, DIST_MAX);
  using P = pair<EdgeTy, NodeId>;
  priority_queue<P, vector<P>, greater<P>> pq;
  dist[s] = 0;
  pq.push({dist[s], s});
  int itr = 0;
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();
    if (u == t) {
      break;
    }
    if (d != dist[u]) {
      continue;
    }
    itr++;
    for (size_t i = G.offset[u]; i < G.offset[u + 1]; i++) {
      NodeId v = G.E[i].v;
      EdgeTy w = G.E[i].w;
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push({dist[v], v});
      }
    }
  }
  return make_pair(dist[t], itr);
}

pair<EdgeTy, int> PchQuery::stQuery(NodeId s, NodeId t) {
  if (s == t) {
    return {0, 0};
  }
  s = GC.rank[s], t = GC.rank[t];
  current_timestamp++;
  using P = pair<EdgeTy, NodeId>;
  priority_queue<P, vector<P>, greater<P>> pq;
  NodeId s_d = s << 1, t_d = t << 1 | 1;
  dist[s_d] = dist[t_d] = 0;
  timestamp[s_d] = timestamp[t_d] = current_timestamp;
  pq.push({dist[s_d], s_d});
  pq.push({dist[t_d], t_d});
  EdgeTy ans = DIST_MAX;
  int itr = 0;

  while (!pq.empty()) {
    auto [d, u_d] = pq.top();
    pq.pop();
    if (d != dist[u_d]) {
      continue;
    }
    if (d >= ans) {
      break;
    }
    itr++;
    NodeId u = u_d >> 1, dir = u_d & 1;
    if (GC.level[u] != 0) {
      const auto &in_offset = dir ? GC.offset : GC.in_offset;
      const auto &in_E = dir ? GC.E : GC.in_E;
      bool stall = false;
      for (EdgeId i = in_offset[u]; i < in_offset[u + 1]; i++) {
        NodeId v_d = in_E[i].v << 1 | dir;
        NodeId w = in_E[i].w;
        if (timestamp[v_d] == current_timestamp) {
          if (dist[v_d] + w <= dist[u_d]) {
            stall = true;
            break;
          }
        }
      }
      if (stall) {
        continue;
      }
    }
    const auto &out_offset = dir ? GC.in_offset : GC.offset;
    const auto &out_E = dir ? GC.in_E : GC.E;
    for (EdgeId i = out_offset[u]; i < out_offset[u + 1]; i++) {
      NodeId v_d = out_E[i].v << 1 | dir;
      NodeId w = out_E[i].w;
      if (timestamp[v_d] != current_timestamp) {
        timestamp[v_d] = current_timestamp;
        dist[v_d] = DIST_MAX;
      }
      if (dist[v_d] > d + w) {
        dist[v_d] = d + w;
        if (dist[v_d] >= ans) {
          continue;
        }
        if (timestamp[v_d ^ 1] == current_timestamp) {
          ans = min(ans, dist[v_d] + dist[v_d ^ 1]);
        }
        // TODO: terminate at dist[v_d]*2 >= ans for uncontractable graph
        pq.push({dist[v_d], v_d});
      }
    }
  }
  return make_pair(ans, itr);
}

// sequence<EdgeTy> PchQuery::ssspQuery(NodeId s) {
//   s = GC.rank[s];
//   current_timestamp++;
//   using P = pair<EdgeTy, NodeId>;
//   priority_queue<P, vector<P>, greater<P>> pq;
//   timestamp[s] = current_timestamp;
//   dist[s] = 0;
//   pq.push({dist[s], s});
//   while (!pq.empty()) {
//     auto [d, u] = pq.top();
//     pq.pop();
//     if (d != dist[u]) {
//       continue;
//     }
//     for (size_t i = GC.offset[u]; i < GC.offset[u + 1]; i++) {
//       NodeId v = GC.E[i].v;
//       EdgeTy w = GC.E[i].w;
//       if (timestamp[v] != current_timestamp) {
//         timestamp[v] = current_timestamp;
//         dist[v] = DIST_MAX;
//       }
//       if (dist[v] > dist[u] + w) {
//         dist[v] = dist[u] + w;
//         pq.push({dist[v], v});
//       }
//     }
//     for (size_t i = in_offset_incoming[u]; i < in_offset_incoming[u + 1];
//     i++) {
//       NodeId v = in_E_incoming[i].v;
//       EdgeTy w = in_E_incoming[i].w;
//       if (timestamp[v] != current_timestamp) {
//         timestamp[v] = current_timestamp;
//         dist[v] = DIST_MAX;
//       }
//       if (dist[v] > dist[u] + w) {
//         dist[v] = dist[u] + w;
//         pq.push({dist[v], v});
//       }
//     }
//   }
//   // TODO: this part could be discounted from the timer
//   sequence<EdgeTy> mapping_dist(n);
//   parallel_for(0, n, [&](NodeId i) {
//     NodeId v = GC.rank[i];
//     if (timestamp[v] != current_timestamp) {
//       mapping_dist[i] = DIST_MAX;
//     } else {
//       mapping_dist[i] = dist[v];
//     }
//   });
//   return mapping_dist;
// }

sequence<EdgeTy> PchQuery::ssspQuery(NodeId s, bool remap = true,
                                     bool in_parallel = false) {
  s = GC.rank[s];
  current_timestamp++;
  using P = pair<EdgeTy, NodeId>;
  priority_queue<P, vector<P>, greater<P>> pq;
  timestamp[s] = current_timestamp;
  dist[s] = 0;
  pq.push({dist[s], s});
  // forward update
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();
    if (d != dist[u]) {
      continue;
    }
    for (size_t i = GC.offset[u]; i < GC.offset[u + 1]; i++) {
      NodeId v = GC.E[i].v;
      EdgeTy w = GC.E[i].w;
      if (timestamp[v] != current_timestamp) {
        timestamp[v] = current_timestamp;
        dist[v] = DIST_MAX;
      }
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push({dist[v], v});
      }
    }
  }
  // backward update
  // int lower_b = 0, upper_b = GC.layerOffset.size() - 1;
  int lower_b = GC.ccOffset[GC.ccRank[s]],
      upper_b = GC.ccOffset[GC.ccRank[s] + 1];
  if (GC.level[GC.layerOffset[lower_b]] == 0) {
    lower_b++;
  }
  if (!in_parallel) {
    for (int i = upper_b - 1; i >= lower_b; --i) {
      // parallel_for(G.layerOffset[i - 1], G.layerOffset[i], [&](size_t k) {
      for (size_t k = GC.layerOffset[i]; k < GC.layerOffset[i + 1]; ++k) {
        if (timestamp[k] != current_timestamp) {
          timestamp[k] = current_timestamp;
          dist[k] = DIST_MAX;
        }
        for (size_t j = GC.in_offset[k]; j < GC.in_offset[k + 1]; j++) {
          NodeId v = GC.in_E[j].v;
          EdgeTy w = GC.in_E[j].w;
          if (timestamp[v] != current_timestamp || dist[v] == DIST_MAX) {
            continue;
          }
          if (dist[k] > dist[v] + w) {
            dist[k] = dist[v] + w;
          }
        }
      }
    }
  } else {
    for (int i = upper_b - 1; i >= lower_b; --i) {
      parallel_for(GC.layerOffset[i], GC.layerOffset[i + 1], [&](size_t k) {
        if (timestamp[k] != current_timestamp) {
          timestamp[k] = current_timestamp;
          dist[k] = DIST_MAX;
        }
        for (size_t j = GC.in_offset[k]; j < GC.in_offset[k + 1]; j++) {
          NodeId v = GC.in_E[j].v;
          EdgeTy w = GC.in_E[j].w;
          if (timestamp[v] != current_timestamp || dist[v] == DIST_MAX) {
            continue;
          }
          if (dist[k] > dist[v] + w) {
            dist[k] = dist[v] + w;
          }
        }
      });
    }
  }
  if (remap) {
    // TODO: this part could be discounted from the timer
    sequence<EdgeTy> mapping_dist(n);
    parallel_for(0, n, [&](NodeId i) {
      NodeId v = GC.rank[i];
      if (timestamp[v] != current_timestamp) {
        mapping_dist[i] = DIST_MAX;
      } else {
        mapping_dist[i] = dist[v];
      }
    });
    return mapping_dist;
  }
  return dist;
}