#include <queue>

#include "graph.hpp"

void dijkstra(size_t s, const Graph &G, EdgeTy *dist) {
  fill(dist, dist + G.n, INT_MAX);
  dist[s] = 0;
  using T = pair<EdgeTy, NodeId>;
  priority_queue<T, vector<T>, greater<T>> pq;
  pq.push(make_pair(dist[s], s));
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();
    if (dist[u] < d) continue;
    for (size_t j = G.offset[u]; j < G.offset[u + 1]; j++) {
      NodeId v = G.E[j].v;
      EdgeTy w = G.E[j].w;
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push(make_pair(dist[v], v));
      }
    }
  }
  return;
}

void decompress_ch(size_t s, const PchGraph &G, EdgeTy *dist,
                   bool in_parallel = false) {
  size_t lower_b = G.ccOffset[G.ccRank[s]],
         upper_b = G.ccOffset[G.ccRank[s] + 1];
  if (G.level[G.layerOffset[lower_b]] == 0) lower_b++;
  if (!in_parallel) {
    for (size_t i = upper_b; i > lower_b; --i) {
      // parallel_for(G.layerOffset[i - 1], G.layerOffset[i], [&](size_t k) {
      for (size_t k = G.layerOffset[i - 1]; k < G.layerOffset[i]; ++k) {
        for (size_t j = G.in_offset[k]; j < G.in_offset[k + 1]; j++) {
          NodeId v = G.in_E[j].v;
          EdgeTy w = G.in_E[j].w;
          if (dist[k] > dist[v] + w) {
            dist[k] = dist[v] + w;
          }
        }
      }  //);
    }
  } else {
    for (size_t i = upper_b; i > lower_b; --i) {
      parallel_for(G.layerOffset[i - 1], G.layerOffset[i], [&](size_t k) {
        for (size_t j = G.in_offset[k]; j < G.in_offset[k + 1]; j++) {
          NodeId v = G.in_E[j].v;
          EdgeTy w = G.in_E[j].w;
          if (dist[k] > dist[v] + w) {
            dist[k] = dist[v] + w;
          }
        }
      });
    }
  }
}

void sssp_correctness_verifier(EdgeTy *cor_dist, EdgeTy *ch_dist,
                               const PchGraph &GC, bool test = false) {
  if (!test) {
    parallel_for(0, GC.n,
                 [&](size_t i) { assert(cor_dist[i] == ch_dist[GC.rank[i]]); });
  } else {
    int itr = 0;
    // size_t layer = 0;
    for (size_t i = 0; i < GC.n; ++i) {
      if (cor_dist[i] != ch_dist[GC.rank[i]]) {
        itr++;
        cout << "cor_dist[" << i << "]=" << cor_dist[i] << ", my_dist[" << i
             << "]=" << ch_dist[GC.rank[i]] << endl;
        if (itr == 300) {
          assert(cor_dist[i] == ch_dist[GC.rank[i]]);
        }
      }
    }
  }
}

void sssp_verifier(const Graph &G, const PchGraph &GC, int source_num = 1,
                   bool verify = true) {
  EdgeTy *cor_dist = new EdgeTy[G.n];
  EdgeTy *ch_dist = new EdgeTy[G.n];
  internal::timer t1;
  internal::timer t2;
  internal::timer t3;

  for (int v = 0; v < source_num; ++v) {
    t1.reset();
    t2.reset();
    t3.reset();
    NodeId s = hash32(v) % G.n;
    t1.start();
    dijkstra(s, G, cor_dist);
    t1.stop();
    t2.start();
    NodeId s1 = GC.rank[s];
    dijkstra(s1, GC, ch_dist);
    sssp_correctness_verifier(cor_dist, ch_dist, GC, true);
    t3.start();
    decompress_ch(s1, GC, ch_dist);
    t2.stop();
    t3.stop();
    if (verify) sssp_correctness_verifier(cor_dist, ch_dist, GC, true);
    printf("%d\t%f\t%f\t%f\n", s, t1.total_time(), t2.total_time(),
           t3.total_time());
  }
  delete[] cor_dist;
  delete[] ch_dist;
}

void test(const Graph &G, const PchGraph &GC, int source_num = 1) {
  EdgeTy *cor_dist = new EdgeTy[G.n];
  EdgeTy *ch_dist = new EdgeTy[G.n];
  internal::timer t1;
  internal::timer t2;
  for (int v = 0; v < source_num; ++v) {
    t1.reset();
    t2.reset();
    NodeId s = hash32(v) % G.n;
    t1.start();
    dijkstra(s, G, cor_dist);
    t1.stop();
    t2.start();
    dijkstra(s, GC, ch_dist);
    t2.stop();
    parallel_for(0, GC.n, [&](size_t i) { assert(cor_dist[i] == ch_dist[i]); });
    printf("%d\t%f\t%f\n", s, t1.total_time(), t2.total_time());
  }
  delete[] cor_dist;
  delete[] ch_dist;
}