#include "graph.hpp"

struct PchQuery {
  size_t n;
  const PchGraph &GC;
  const Graph &G;
  sequence<EdgeTy> dist;
  sequence<uint32_t> timestamp;
  uint32_t current_timestamp;

  PchQuery(const PchGraph &_GC, const Graph &_G) : GC(_GC), G(_G) {
    assert(G.n == GC.n);
    n = G.n;
    dist = sequence<EdgeTy>(n * 2);
    timestamp = sequence<uint32_t>(n * 2);
    current_timestamp = 0;
  }

  EdgeTy stVerifier(NodeId s, NodeId t);
  pair<EdgeTy, int> stQuery(NodeId s, NodeId t);
};

EdgeTy PchQuery::stVerifier(NodeId s, NodeId t) {
  // TODO: change to bidirectional search
  sequence<EdgeTy> exp_dist(n, DIST_MAX);
  using P = pair<EdgeTy, NodeId>;
  priority_queue<P, vector<P>, greater<P>> pq;
  exp_dist[s] = 0;
  pq.push({exp_dist[s], s});
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();
    if (u == t) {
      break;
    }
    if (d != exp_dist[u]) {
      continue;
    }
    for (size_t i = G.offset[u]; i < G.offset[u + 1]; i++) {
      NodeId v = G.E[i].v;
      EdgeTy w = G.E[i].w;
      if (exp_dist[v] > exp_dist[u] + w) {
        exp_dist[v] = exp_dist[u] + w;
        pq.push({exp_dist[v], v});
      }
    }
  }
  return exp_dist[t];
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
    itr++;
    auto [d, u_d] = pq.top();
    pq.pop();
    if (d != dist[u_d]) {
      continue;
    }
    if (d >= ans) {
      break;
    }
    NodeId u = u_d >> 1, dir = u_d & 1;
    const auto &out_offset = dir ? GC.in_offset : GC.offset;
    const auto &out_E = dir ? GC.in_E : GC.E;
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
