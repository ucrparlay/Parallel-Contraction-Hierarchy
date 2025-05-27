SPoCH: Scalable Parallelization of Contraction Hierarchies
====================== 

This repository includes the implmentation of SPoCH (Scalable Parallelization of Contraction Hierarchies). 

## Developing 

### Prerequisites 
* g++ or clang with C++17 features support (tested with clang 16.0.6 and g++ 14.2.1) on Linux machines.
* We use [ParlayLib](https://github.com/cmuparlay/parlaylib) to support fork-join parallelism and some parallel primitives. It is provided as a submodule in our repository. 
### Setting up 
Clone the library with submodule 
```shell
git clone --recurse-submodules https://github.com/ucrparlay/Parallel-Contraction-Hierarchy.git
cd Parallel-Contraction-Hierarchy/ 
```

Alternatively, you can first clone it and add the submodule 
```shell
git clone https://github.com/ucrparlay/Parallel-Contraction-Hierarchy.git
git submodule update --init --recursive 
cd Parallel-Contraction-Hierarchy/ 
```

### Building
A makefile is given in the repository, you can compile the code by: 
```shell
make
```

## Usage

### Build Parallel Contraction Hierarchy
```shell
./build_pch [-i input_file_path] [-o output_file_path] [-p max_pop_count] [-s selection_fraction] [-t st_query_verify_num] [-q sssp_query_verify_num] [-d print_details] 
```
Options: 
* `-i <input_file_path>`: (Required) the path to the input file containing graph data. Feasible graphs can be found [here.](https://pasgal-bs.cs.ucr.edu/)
* `-o <output_file_path>`: the path to save the output results.
* `-p <max_pop_count>`: the maximum number of vertices to settle during the local search phase. Higher values may improve accuracy at the cost of runtime.
* `-s <selection_fraction>`: the fraction of feasible vertices to be contracted in each round. For example, 0.1 selects 10% of feasible vertices.
* `-t <st_query_verify_num>`: the number of s-t queries to verify the correctness of output CH graph.
* `-q <sssp_query_verify_num>`: the number of sssp queries to verify the correctness of output CH graph.
* `-d` print the details of each round

For example, if you want to compute the CH of a weighted graph INPUT_NAME and store the output CH graph in OUTPUT_NAME, set max_pop_count=5000, selection_fraction=1, run st query 20 times and print the per round detail during the contraction, you can run:
```shell
./build_pch -i INPUT_NAME -o OUTPUT_NAME -p 5000 -s 1 -t 20 -d
```
## Graph Formats
The application can auto-detect the format of the input graph based on the suffix of the filename. Here is a list of supported graph formats: 
+ `.bin` The binary graph format from [GBBS](https://github.com/ParAlg/gbbs). 
+ `.adj` The adjacency graph format from [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). 


### Query

```shell
./query [origin_graph] [ch_graph] [st_query_verify_num] [sssp_query_verify_num]
```

For example, for a origin_graph ORIGIN_GRAPH and its corresponding ch_graph CH_GRAPH, if you want to run st query 20 times and sssp query 40 times, you can run:
```shell
./query ORIGIN_GRAPH CH_GRAPH 20 40
```

Reference
--------
Zijin Wan, Xiaojun Dong, Letong Wang, Enzuo zhu, Yan Gu and Yihan Sun. [*Parallel Contraction Hierarchies Can Be Efficient and Scalable*](https://arxiv.org/abs/2412.18008). To appear in International Conference on Supercomputing (ICS), 2025

If you use our code, please cite our paper:
```
@inproceedings{wan2025parallel,
  author    = {Wan, Zijin and Dong, Xiaojun and Wang, Letong and Zhu, Enzuo and Gu, Yan and Sun, Yihan},
  title     = {Parallel Contraction Hierarchies Can Be Efficient and Scalable},
  booktitle = {International Conference on Supercomputing (ICS)},
  year      = {2025}
}
```
