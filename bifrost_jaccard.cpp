#include <algorithm>
#include <bifrost/ColoredCDBG.hpp>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

int main(int argc, char **argv) {
  const std::string filename_prefix(argv[1]);
  const std::string colors_filename(argv[2]);
  const std::string table_filename(argv[3]);

  // open graph with colors
  ColoredCDBG<> cdbg;
  cdbg.read(filename_prefix + ".gfa.gz", filename_prefix + ".bfi",
            filename_prefix + ".color.bfg");

  // output color list
  std::ofstream cofs(colors_filename);
  for (const auto &c : cdbg.getColorNames())
    cofs << c << std::endl;
  cofs.close();

  // store k-mer lists for each color
  std::unordered_map<size_t, std::unordered_set<size_t>> ckmap;

  std::hash<string> hasher;
  int i = 0;
  int k = 0;
  for (const UnitigColorMap<void> &unitig : cdbg) {
    UnitigColors *colors = unitig.getData()->getUnitigColors(unitig);
    for (auto it = colors->begin(unitig); it != colors->end(); it++) {
      size_t kmer_pos = it.getKmerPosition();
      size_t color_id = it.getColorID();

      // create k-mer set if color not in map
      if (ckmap.find(color_id) == ckmap.end()) {
        ckmap.emplace(std::make_pair(color_id, std::unordered_set<size_t>()));
      }
      ckmap[color_id].insert(hasher(unitig.getUnitigKmer(kmer_pos).toString()));

      k = k + 1;
    }
    if (i % 100 == 0)
      std::cerr << "Colors: " << ckmap.size() << " , k-mers: " << k
                << std::endl;
    i = i + 1;
  }

  std::ofstream table_ofs(table_filename);

  // header
  table_ofs << "color1" << "\t" << "color2" << "\t" <<
    "color1_kmers" << "\t" << "color2_kmers" << "\t" <<
    "intersect" << "\t" << "union" << "\t" <<
    "jaccard" << std::endl;

  // pairwise jaccard distances
  for (auto c1 : ckmap) {
    for (auto c2 : ckmap) {
      std::vector<size_t> color_intersect;
      std::unordered_set<size_t> color_union;

      // intersection
      for (const auto &k : c1.second) {
        if (c2.second.find(k) != c2.second.end()) {
          color_intersect.push_back(k);
        }
      }

      // union
      for (const auto &k : c1.second) {
        color_union.insert(k);
      }
      for (const auto &k : c2.second) {
        color_union.insert(k);
      }

      table_ofs << c1.first << "\t" << c2.first << "\t" << c1.second.size()
                << "\t" << c2.second.size() << "\t" << color_intersect.size()
                << "\t" << color_union.size() << "\t"
                << (float)color_intersect.size() / color_union.size()
                << std::endl;
    }
  }
  table_ofs.close();
}
