#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

// Containers for the local searches of PCH::pruneNeighbors. A search runs for
// every overlay vertex in a round; with std::unordered_map / unordered_set
// (one heap allocation per insertion) and a fresh priority_queue per search,
// allocation dominated the running time. These containers keep their memory:
// each thread owns one instance, which clear() resets in O(1) for the next
// search.

// Open addressing with linear probing. clear() does not touch the entries: it
// starts a new round, and an entry is present only if it was written in the
// current round.
template <typename K, typename V>
struct flat_hash_map {
  std::vector<K> keys;
  std::vector<V> values;
  std::vector<uint32_t> rounds;
  uint32_t round = 0;
  size_t num_keys = 0;

  void clear() {
    if (keys.empty()) {
      reset(1 << 10);
    } else if (++round == 0) {  // the round counter wrapped around
      std::fill(rounds.begin(), rounds.end(), 0);
      round = 1;
    }
    num_keys = 0;
  }
  bool contains(K k) const { return rounds[slot(k)] == round; }
  const V *find(K k) const {  // nullptr if k is absent
    size_t i = slot(k);
    return rounds[i] == round ? &values[i] : nullptr;
  }
  V *find(K k) {
    size_t i = slot(k);
    return rounds[i] == round ? &values[i] : nullptr;
  }
  V &operator[](K k) {  // inserts V() if k is absent, as std::unordered_map
    size_t i = slot(k);
    if (rounds[i] == round) return values[i];
    if (2 * (num_keys + 1) > keys.size()) {  // keep the load factor <= 1/2
      grow();
      i = slot(k);
    }
    keys[i] = k;
    values[i] = V();
    rounds[i] = round;
    num_keys++;
    return values[i];
  }

 private:
  // the entry of k, or the free entry where k would be inserted
  size_t slot(K k) const {
    size_t mask = keys.size() - 1;
    size_t i = (static_cast<size_t>(k) * 0x9E3779B97F4A7C15ull >> 20) & mask;
    while (rounds[i] == round && keys[i] != k) i = (i + 1) & mask;
    return i;
  }
  void reset(size_t capacity) {  // capacity: a power of two
    keys.assign(capacity, K());
    values.assign(capacity, V());
    rounds.assign(capacity, 0);
    round = 1;
    num_keys = 0;
  }
  void grow() {  // doubles the capacity and reinserts the present entries
    std::vector<K> old_keys = std::move(keys);
    std::vector<V> old_values = std::move(values);
    std::vector<uint32_t> old_rounds = std::move(rounds);
    uint32_t old_round = round;
    reset(2 * old_keys.size());
    for (size_t i = 0; i < old_keys.size(); i++) {
      if (old_rounds[i] != old_round) continue;
      size_t j = slot(old_keys[i]);
      keys[j] = old_keys[i];
      values[j] = old_values[i];
      rounds[j] = round;
      num_keys++;
    }
  }
};

// A set on top of flat_hash_map; erase() only marks the entry as absent.
template <typename K>
struct flat_hash_set {
  flat_hash_map<K, uint8_t> map;  // 1: k is in the set (0: erased)
  size_t size = 0;

  void clear() {
    map.clear();
    size = 0;
  }
  bool contains(K k) const {
    const uint8_t *in = map.find(k);
    return in && *in;
  }
  bool empty() const { return size == 0; }
  void insert(K k) {
    uint8_t &in = map[k];
    if (!in) {
      in = 1;
      size++;
    }
  }
  void erase(K k) {
    uint8_t *in = map.find(k);
    if (in && *in) {
      *in = 0;
      size--;
    }
  }
};

// std::priority_queue (smallest first) whose storage is kept by clear(), so
// the pop order is exactly that of std::priority_queue.
template <typename T>
struct min_heap : std::priority_queue<T, std::vector<T>, std::greater<T>> {
  void clear() { this->c.clear(); }
};
