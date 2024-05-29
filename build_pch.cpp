#include "build_pch.hpp"
// #include "dijkstra.hpp"
#include "query.hpp"

int main(int arcontracted_graph, char *argv[]) {
  if (arcontracted_graph < 3) {
    fprintf(stderr,
            "WeightedSymmetricGraph only\n"
            "Usage: %s [input_file] [output_file]\nOptional: [max_pop_count] "
            "[bidirect_verify_num] [sssp_verify_num]\n",
            argv[0]);
    return 0;
  }
  INPUT_FILEPATH = argv[1];
  OUTPUT_FILEPATH = argv[2];
  size_t max_pop_count = 500;
  int bidirect_verify_num = 0;
  int sssp_verify_num = 0;
  if (arcontracted_graph > 3) max_pop_count = atol(argv[3]);
  if (arcontracted_graph > 4) bidirect_verify_num = atol(argv[4]);
  if (arcontracted_graph > 5) sssp_verify_num = atol(argv[5]);

  // #ifdef DEBUG
  //   size_t n = 100, m = 260;
  //   for (int seed = 0; seed <= 1000000; seed++) {
  //     printf("seed: %d\n", seed);
  //     Graph origin_graph = generate_random_graph(n, m, seed);
  //     PCH *solver = new PCH(origin_graph, max_pop_count);
  //     PchGraph contracted_graph = solver->createContractionHierarchy();
  //     delete (solver);
  //     PchQuery query(contracted_graph, origin_graph);
  //     ofstream ofs("pch.tsv");
  //     ofs << fixed << setprecision(6);
  //     printf("Start query\n");
  //     for (int i = 0; i < bidirect_verify_num; i++) {
  //       NodeId s = hash32(i) % origin_graph.n;
  //       NodeId t = hash32(i + s) % origin_graph.n;
  //       internal::timer tm;
  //       EdgeTy d = query.stQuery(s, t);
  //       tm.stop();
  //       // if (d == DIST_MAX) {
  //       //   continue;
  //       // }
  //       EdgeTy exp_dist = query.stVerifier(s, t);
  //       if (exp_dist != d) {
  //         printf("Error: s: %u, t: %u, output_dist: %u, exp_dist: %u\n", s,
  //         t, d,
  //                exp_dist);
  //         abort();
  //       }
  //       // printf("s: %u, t: %u, d: %u\n", s, t, d);
  //       ofs << s << '\t' << t << '\t' << d << '\t' << tm.total_time() <<
  //       '\n';
  //     }
  //     ofs.close();
  //   }
  //   return 0;
  // #endif

  Graph origin_graph = read_graph(INPUT_FILEPATH);
  //   for (size_t i = 0; i < origin_graph.n; i++) {
  //     printf("N(%zu): ", i);
  //     for (size_t j = origin_graph.offset[i]; j < origin_graph.offset[i + 1];
  //          j++) {
  //       printf("(%d,%d) ", origin_graph.E[j].v, origin_graph.E[j].w);
  //     }
  //     puts("");
  //   }
  PCH *solver = new PCH(origin_graph, max_pop_count);
  PchGraph contracted_graph = solver->createContractionHierarchy();
  delete (solver);
  PchQuery query(contracted_graph, origin_graph);
  ofstream ofs("pch.tsv");
  ofs << fixed << setprecision(6);
  printf("Start query\n");
  for (int i = 0; i < bidirect_verify_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    NodeId t = hash32(i + s) % origin_graph.n;
    internal::timer tm;
    auto [d, itr] = query.stQuery(s, t);
    tm.stop();
    // if (d == DIST_MAX) {
    //   continue;
    // }
    // EdgeTy exp_dist = query.stVerifier(s, t);
    // if (exp_dist != d) {
    //   printf("Error: s: %u, t: %u, output_dist: %u, exp_dist: %u\n", s, t, d,
    //          exp_dist);
    //   abort();
    // }
    // printf("s: %u, t: %u, d: %u\n", s, t, d);
    ofs << s << '\t' << t << '\t' << itr << '\t' << d << '\t' << tm.total_time()
        << '\n';
  }
  ofs.close();
  //   sssp_verifier(origin_graph, contracted_graph, bidirect_verify_num);
  //   sssp_verifier(origin_graph, contracted_graph, sssp_verify_num);
  return 0;
}
