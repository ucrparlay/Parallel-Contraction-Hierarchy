#pragma once
#include "graph.hpp"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "utilities.h"

inline NodeId find_compress(NodeId i, parlay::sequence<NodeId> &parents) {
  NodeId j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  NodeId tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

inline bool unite_impl(NodeId u_orig, NodeId v_orig,
                       parlay::sequence<NodeId> &parents) {
  NodeId u = u_orig;
  NodeId v = v_orig;
  while (u != v) {
    u = find_compress(u, parents);
    v = find_compress(v, parents);
    if (u > v && parents[u] == u &&
        atomic_compare_and_swap(&parents[u], u, v)) {
      return true;
    } else if (v > u && parents[v] == v &&
               atomic_compare_and_swap(&parents[v], v, u)) {
      return true;
    }
  }
  return false;
}

// Returns the connected component id of each vertex
template <class Graph>
parlay::sequence<NodeId> get_cc(const Graph &G) {
  size_t n = G.n;
  parlay::sequence<NodeId> parents(n);
  parlay::parallel_for(0, n, [&](size_t i) { parents[i] = i; });
  parlay::parallel_for(0, n, [&](NodeId u) {
    parlay::parallel_for(
        G.offset[u], G.offset[u + 1],
        [&](size_t i) {
          NodeId v = G.E[i].v;
          unite_impl(u, v, parents);
        },
        1024);
  });
  parlay::parallel_for(
      0, n, [&](NodeId i) { parents[i] = find_compress(i, parents); });
  return parents;
}
