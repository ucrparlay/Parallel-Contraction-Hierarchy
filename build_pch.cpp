#include "build_pch.hpp"

#include "dijkstra.hpp"
#include "query.hpp"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr,
            "WeightedGraph only\n"
            "Usage: %s [-i input_graph]\n"
            "Options:\n"
            "\t-o,\toutput_ch_graph\n"
            "\t-p,\tmax_pop_count\n"
            "\t-s,\tselect fraction\n"
            "\t-t,\ts-t query verify numn"
            "\t-q,\tsssp query verify num\n"
            "\t-b,\tdegree bound\n"
            "\t-d,\tprint per round detail\n",
            argv[0]);
    return 0;
  }
  char c;
  size_t max_pop_count = 500;
  int bidirect_verify_num = 0;
  int sssp_verify_num = 0;
  NodeId degree_bound = 0;
  bool degree_bounded = false, print_detail = false, write_ch = false;
  double sample_bound = 1;
  while ((c = getopt(argc, argv, "i:o:p:s:t:q:b:d")) != -1) {
    switch (c) {
      case 'i':
        INPUT_FILEPATH = optarg;
        break;
      case 'o':
        write_ch= true;
        OUTPUT_FILEPATH = optarg;
        break;
      case 'p':
        max_pop_count = atol(optarg);
        if (max_pop_count < 0) {
          fprintf(stderr, "Error: max_pop_count must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 's':
        sample_bound = atof(optarg);
        if (sample_bound <= 0 || sample_bound >1) {
          fprintf(stderr, "Error: selection_fraction must be larger than 0 and smaller than or equal to 1\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 't':
        bidirect_verify_num = atol(optarg);
        if (bidirect_verify_num < 0) {
          fprintf(stderr, "Error: bidirect_verify_num must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'q':
        sssp_verify_num = atol(optarg);
        if (sssp_verify_num < 0) {
          fprintf(stderr, "Error: sssp_verify_num must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'b':
        degree_bounded = true;
        degree_bound = atol(optarg);
        if (degree_bound < 0) {
          fprintf(stderr, "Error: degree_bound must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'd':
        print_detail = true;
        break;
      default:
        fprintf(stderr, "Error: Unknown option %c\n", optopt);
        exit(EXIT_FAILURE);
    }
  }
  Graph origin_graph = read_graph(INPUT_FILEPATH);
  PCH *solver =
      new PCH(origin_graph, max_pop_count, degree_bounded, degree_bound, sample_bound, print_detail);
  PchGraph contracted_graph = solver->createContractionHierarchy();
  delete (solver);
  PchQuery query(contracted_graph, origin_graph);
  ofstream ofs("pch.tsv", ios::app);
  ofs << fixed << setprecision(6);
  printf("Start query\n");
  query.make_inverse();
  for (int i = 0; i < bidirect_verify_num; i++) {
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
  for (int i = 0; i < sssp_verify_num; i++) {
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
  ofs.close();
  if(write_ch)write_pbbs_format(contracted_graph, OUTPUT_FILEPATH);
  return 0;
}
