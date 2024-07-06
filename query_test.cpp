#include "query.hpp"

#include "graph.hpp"
#include "iostream"

int main(int argc, char *argv[]) {
  if (argc < 5) {
    fprintf(stderr,
            "WeightedSymmetricGraph only\n"
            "Usage: %s [origin_graph] [contracted_graph] [bidirect_query_num] "
            "[bidirect_query_num]\n",
            argv[0]);
    return 0;
  }
  char *FILEPATH = argv[1];
  char *FILEPATH2 = argv[2];
  int bidirect_query_num = atol(argv[3]);
  int sssp_query_num = atol(argv[4]);
  Graph origin_graph = read_graph(FILEPATH);
  PchGraph contracted_graph;
  read_contracted_graph(FILEPATH2, contracted_graph);
  cout << "finish" << endl;
  PchQuery query(contracted_graph, origin_graph);
  ofstream ofs("query.tsv", ios::app);
  ofs << fixed << setprecision(6);
  for (int i = 0; i < bidirect_query_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    NodeId t = hash32(i + s) % origin_graph.n;
    internal::timer tm;
    auto [d, itr] = query.stQuery(s, t);
    tm.stop();
    internal::timer tm2;
    auto [exp_dist, itr2] = query.stVerifier(s, t);
    tm2.stop();
    if (exp_dist != d) {
      printf("Error: s: %u, t: %u, output_dist: %u, exp_dist: %u\n", s, t, d,
             exp_dist);
      abort();
    }
    printf("s: %u, t: %u, itr: %u, d: %u\n", s, t, itr, d);
    ofs << s << '\t' << t << '\t' << itr << '\t' << itr2 << '\t' << d << '\t'
        << tm.total_time() << '\t' << tm2.total_time() << '\n';
  }
  for (int i = 0; i < sssp_query_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    internal::timer tm;
    auto dist = query.ssspQuery(s);
    tm.stop();
    printf("%d SSSP query: s: %u\n", i, s);
    internal::timer tm2;
    auto exp_dist = query.ssspVerifier(s);
    tm2.stop();
    if (exp_dist != dist) {
      for (size_t i = 0; i < dist.size(); i++) {
        if (dist[i] != exp_dist[i]) {
          printf("output_dist[%zu]: %u, exp_dist[%zu]: %u\n", i, dist[i], i,
                 exp_dist[i]);
        }
      }
      abort();
    }
    ofs << s << '\t' << tm.total_time() << '\t' << tm2.total_time() << '\n';
  }
  // for (int i = 0; i < bidirect_query_num; i++) {
  //   NodeId s = hash32(i) % origin_graph.n;
  //   NodeId t = hash32(i + s) % origin_graph.n;
  //   internal::timer tm;
  //   auto [d, itr] = query.stQuery(s, t);
  //   tm.stop();
  //   printf("s: %u t: %u itr: %u d: %u, time: %f\n", s, t, itr, d,
  //          tm.total_time());
  //   ofs << s << '\t' << t << '\t' << itr << '\t' << d << '\t' << tm.total_time()
  //       << '\n';
  // }
  // for (int i = 0; i < sssp_query_num; i++) {
  //   NodeId s = hash32(i) % origin_graph.n;
  //   internal::timer tm;
  //   auto dist = query.ssspQuery(s, false, true);
  //   tm.stop();
  //   printf("s: %u, time: %f\n", s, tm.total_time());
  //   ofs << s << '\t' << tm.total_time() << '\n';
  // }
  ofs.close();
  return 0;
}
