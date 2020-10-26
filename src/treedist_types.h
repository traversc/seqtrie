#include <Rcpp.h>

struct SeqNode;
using SeqNodeUptr = std::unique_ptr<SeqNode>;

struct SeqNode {
  std::set<int> idx;
  std::unordered_map<char, SeqNodeUptr > seqmap;
  SeqNode() : idx(std::set<int>()), seqmap(std::unordered_map<char, SeqNodeUptr >()) {}
  bool isLeaf() const {
    if(seqmap.size() > 0) {
      return false;
    } else {
      return true;
    }
  }
  std::set<int> allChildIdx() const {
    if(isLeaf()) {
      return idx;
    } else {
      std::set<int> ret = idx;
      for(auto & e : seqmap) {
        std::set<int> x = e.second->allChildIdx();
        ret.insert(x.begin(), x.end());
      }
      return ret;
    }
  }
};
