#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <optional>
#include <utility>

#include "utilities.h"

// neighbor idx, wgh and hop
template <typename K, typename V1, typename V2, typename Hash = parlay::hash<K>,
          typename Equal = std::equal_to<>>
struct hash_map {
 public:
  static constexpr K KEY_MAX = numeric_limits<K>::max();
  static constexpr V1 V1_MAX = numeric_limits<V1>::max();
  static constexpr V2 V2_MAX = numeric_limits<V2>::max();
  using V = std::pair<V1, V2>;
  using KV = std::pair<K, V>;
  using index = K;
  long m;
  index first_index(K k) const { return hash(k + 1) % m; }
  index next_index(index h) const { return (h + 1 == m) ? 0 : h + 1; }
  Hash hash;
  Equal equal;

  enum state : char { empty, full, locked };
  struct entry {
    K key;
    V value;
    V1 hop;
    bool valid;
    bool settled;
    entry()
        : key(KEY_MAX),
          value(V1_MAX, V2_MAX),
          hop(1),
          valid(true),
          settled(false) {}
  };
  parlay::sequence<entry> H;

  hash_map(long size, Hash &&hash = {}, Equal &&equal = {})
      : m(100 + static_cast<long>(1.5 * size)),
        hash(hash),
        equal(equal),
        H(parlay::sequence<entry>(m)) {}

  std::pair<bool, K> insert(const K &k, const V &v, const V1 h) {
    index i = first_index(k);
    assert(k != v.first);
    while (true) {
      if (H[i].valid == true) {
        if (H[i].key == KEY_MAX) {
          // H[i].valid = true;
          CAS(&H[i].key, KEY_MAX, k);
        }
        if (H[i].key == k) {
          bool flag = false;
          if (CAS(&H[i].value.first, V1_MAX, v.first)) {
            flag = true;
          }
          if (H[i].value.first == v.first) {
            // if(H[i].value.first == 646704) {
            //     printf("adding %u to 646704\n", H[i].key);
            // }
            if (write_min(&H[i].value.second, v.second, std::less<V2>())) {
              H[i].hop = h;
            }
            return make_pair(flag, i);
          }
        }
      }
      i = next_index(i);
    }
  }

  std::optional<V> find(const K &k) {
    index i = first_index(k);
    while (true) {
      if (H[i].key == KEY_MAX) return {};
      if (H[i].key == k) return H[i].value;
      i = next_index(i);
    }
  }

  bool delete_index(const K &k) {
    // if(H[k].key==KEY_MAX)return true;
    H[k].valid = true;
    H[k].settled = false;
    H[k].key = KEY_MAX;
    H[k].value.first = V1_MAX;
    H[k].value.second = V2_MAX;
    return true;
  }

  void clear_all() {
    parallel_for(0, m, [&](K i) { delete_index(i); });
  }

  parlay::sequence<K> keys() {
    return parlay::map_maybe(H, [](const entry &x) {
      return (x.status == full) ? std::optional{x.key} : std::optional<K>{};
    });
  }

  size_t size() {
    return parlay::reduce(parlay::delayed_map(H, [&](const entry &x) -> long {
      return x.key != KEY_MAX && x.valid == true;
    }));
  }
};

// neighbor idx and wgh only
template <typename K1, typename K2, typename V,
          typename Hash = parlay::hash<uint64_t>,
          typename Equal = std::equal_to<>>
struct hash_map2 {
 public:
  static constexpr K1 KEY_MAX = numeric_limits<K1>::max();
  static constexpr V V_MAX = numeric_limits<V>::max();
  using K = std::pair<K1, K2>;
  using KV = std::pair<K, V>;
  using index = K1;
  uint64_t m;
  index first_index(K k) const {
    return hash(((uint64_t)k.first << 32 | k.second) + 1) % m;
  }
  index next_index(index h) const { return (h + 1 == m) ? 0 : h + 1; }
  Hash hash;
  Equal equal;

  enum state : char { empty, full, locked };
  struct entry {
    K key;
    V value;
    entry() : key(KEY_MAX, KEY_MAX), value(V_MAX) {}
  };
  parlay::sequence<entry> H;

  hash_map2(long size, Hash &&hash = {}, Equal &&equal = {})
      : m(100 + static_cast<long>(1.5 * size)),
        hash(hash),
        equal(equal),
        H(parlay::sequence<entry>(m)) {}

  void insert(K &k, V v) {
    if (k.first == k.second) return;
    index i = first_index(k);
    while (true) {
      if (H[i].key.first == KEY_MAX) {
        CAS(&H[i].key.first, KEY_MAX, k.first);
      }
      if (H[i].key.first == k.first) {
        if (H[i].key.second == KEY_MAX) {
          CAS(&H[i].key.second, KEY_MAX, k.second);
        }
        if (H[i].key.second == k.second) {
          write_min(&H[i].value, v, std::less<V>());
          return;
        }
      }
      i = next_index(i);
    }
  }

  V find(K &k) {
    if (k.first == k.second) return 0;
    index i = first_index(k);
    while (true) {
      if (H[i].key.first == KEY_MAX) {
        return V_MAX;
      }
      if (H[i].key.first == k.first && H[i].key.second == k.second) {
        return H[i].value;
      }
      i = next_index(i);
    }
  }
};
