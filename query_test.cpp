#include "graph.hpp"
#include "query.hpp"

int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    fprintf(
        stderr,
        "WeightedSymmetricGraph only\n"
        "Usage: %s [origin_graph] [contracted_graph] [bidirect_query_num] [bidirect_query_num]\n",
        argv[0]);
    return 0;
  }
  char *FILEPATH = argv[1];
  char *FILEPATH2 = argv[2];
  int bidirect_query_num = atol(argv[3]);
  int sssp_query_num = atol(argv[4]);
  Graph origin_graph = read_graph(FILEPATH);
  PchGraph contracted_graph;
  read_contracted_graph(FILEPATH2,contracted_graph);
  cerr<<"finish"<<endl;
  PchQuery query(contracted_graph, origin_graph);
  for (int i = 0; i < bidirect_query_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    NodeId t = hash32(i + s) % origin_graph.n;
    internal::timer tm;
    auto [d, itr] = query.stQuery(s, t);
    tm.stop();
    printf("%u\t%u\t%u\n", s, t, d);
  }
  for (int i = 0; i < sssp_query_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    internal::timer tm;
    auto dist = query.ssspQuery(s);
    tm.stop();
    cout << s << '\t' << tm.total_time() << '\n';
  }
  return 0;
}
