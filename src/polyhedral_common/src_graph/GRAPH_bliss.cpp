#include "GRAPH_bliss.h"

bliss::Graph* ReadGraphFromFile(FILE *f, unsigned int &nof_vertices)
{
  unsigned int nof_edges;
  int ret;
  bliss::Graph *g =0;
  ret=fscanf(f, "%u %u\n", &nof_vertices, &nof_edges);
  if (ret != 1) {
    std::cerr << "fscanf error while reading graph 1\n";
    throw TerminalException{1};
  }
  g = new bliss::Graph(nof_vertices);
  for (int i=0; i<int(nof_vertices); i++) {
    unsigned int color;
    ret=fscanf(f, "%u\n", &color);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 2\n";
      throw TerminalException{1};
    }
    g->change_color(i, color);
  }
  for (int iEdge=0; iEdge<int(nof_edges); iEdge++) {
    int a, b;
    ret=fscanf(f, "%u %u\n", &a, &b);
    if (ret != 1) {
      std::cerr << "fscanf error while reading graph 3\n";
      throw TerminalException{1};
    }
    g->add_edge(a-1, b-1);
  }
  return g;
}




