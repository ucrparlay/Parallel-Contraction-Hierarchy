#include "query.hpp"
#include "graph.hpp"
#include "iostream"
int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr,
            "WeightedSymmetricGraph only\n"
            "Usage: %s [origin_graph] [bidirect_query_num] "
            "[bidirect_query_num]\n",
            argv[0]);
    return 0;
  }
  char *FILEPATH = argv[1];
  int bidirect_query_num = atol(argv[2]);
  Graph origin_graph = read_graph(FILEPATH);
  PchGraph contracted_graph;
  contracted_graph.n = origin_graph.n;
  cout << "finish" << endl;
  PchQuery query(contracted_graph, origin_graph);
  ofstream ofs("query.tsv", ios::app);
  ofs << fixed << setprecision(6);
  for (int i = 0; i < bidirect_query_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    NodeId t = hash32(i + s) % origin_graph.n;
    internal::timer tm;
    auto [d, itr] = query.stVerifier(s, t);
    tm.stop();
    printf("s: %u t: %u itr: %u d: %u, time: %f\n", s, t, itr, d,
           tm.total_time());
    ofs << s << '\t' << t << '\t' << itr << '\t' << d << '\t' << tm.total_time()
        << '\n';
  }
  ofs.close();
  return 0;
}
