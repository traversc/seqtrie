#include <set>
#include <memory>
#include <tuple>
#include <limits.h> // INT_MAX

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]

// helper functions
std::vector<int> int_range(int start, int end) {
  if(end < start) {
    throw std::runtime_error("end must be larger or equal start");
  }
  std::vector<int> ret(end - start + 1);
  for(int i=0; i<end - start + 1; ++i) {
    ret[i] = start + i;
  }
  return ret;
}
template <typename T> static void set_insert(std::set<uint32_t> & x, const T & new_x) {
  static_assert(std::is_same_v<uint32_t, T> || std::is_same_v<std::set<uint32_t>, T>);
  if constexpr (std::is_same_v<uint32_t, T>) {
    x.insert(new_x);
  } else if constexpr (std::is_same_v<std::set<uint32_t>, T>) {
    x.insert(new_x.begin(), new_x.end());
  }
}


using ustring = std::basic_string<uint8_t>;
using ustring_view = std::basic_string_view<uint8_t>;

template <uint8_t N> struct RadixNode;
template <uint8_t N> using RadixNodeUptr = std::unique_ptr<RadixNode<N>>;

template <uint8_t N> struct RadixNode {
  ustring branch;
  std::array<RadixNodeUptr<N>, N> child_nodes;
  std::set<uint32_t> terminal_idx;
  RadixNode() {}
  bool isLeaf() const {
    for(size_t i=0; i<N; ++i) {
      if(child_nodes[i]) return false;
    }
    return true;
  }
  size_t size() const { // number of children
    size_t x;
    for(size_t i=0; i<N; ++i) {
      if(child_nodes[i]) ++x;
    }
    return x;
  }
  
  template <typename T> static void insert(RadixNodeUptr<N> & node, const ustring_view sequence, const T & seq_idx) {
    if(sequence.size() == 0) {
      // std::cout << "case -1: end of sequence" << std::endl;
      set_insert<T>(node->terminal_idx, seq_idx);
      return;
    }
    uint8_t s = sequence[0];
    if(!node->child_nodes[s]) {
      // std::cout << "case 0: no existing branch" << std::endl;
      node->child_nodes[s].reset(new RadixNode);
      node->child_nodes[s]->branch = sequence;
      set_insert<T>(node->child_nodes[s]->terminal_idx, seq_idx);
      return;
    }
    size_t i = 0;
    while(i < node->child_nodes[s]->branch.size() && i < sequence.size() && sequence[i] == node->child_nodes[s]->branch[i]) { ++i; }
    
    if(i == sequence.size() && i == node->child_nodes[s]->branch.size()) {
      // std::cout << "case 1: sequence is same as branch" << std::endl;
      set_insert<T>(node->child_nodes[s]->terminal_idx, seq_idx);
    } else if(i == sequence.size()) {
      // std::cout << "case 2: sequence is shorter than branch" << std::endl;
      ustring branch_prefix = node->child_nodes[s]->branch.substr(0,i);
      ustring branch_suffix = node->child_nodes[s]->branch.substr(i);
      uint8_t s_insert = branch_suffix[0];
      
      RadixNodeUptr<N> inserted_node(new RadixNode);
      inserted_node->child_nodes[s_insert] = std::move(node->child_nodes[s]);
      inserted_node->branch = std::move(branch_prefix);
      set_insert<T>(inserted_node->terminal_idx, seq_idx);
      inserted_node->child_nodes[s_insert]->branch = std::move(branch_suffix);
      node->child_nodes[s] = std::move(inserted_node);
      
    } else if(i == node->child_nodes[s]->branch.size()) {
      // std::cout << "case 3: sequence is longer than branch" << std::endl;
      ustring_view seq_suffix = sequence.substr(i); // will remain as ustring_view
      insert<T>(node->child_nodes[s], seq_suffix, seq_idx);
    } else {
      // std::cout << "case 4: sequence is different than branch" << std::endl;
      ustring branch_prefix = node->child_nodes[s]->branch.substr(0,i);
      ustring branch_suffix = node->child_nodes[s]->branch.substr(i);
      ustring seq_suffix = ustring(sequence.substr(i));
      uint8_t s_insert_branch = branch_suffix[0];
      uint8_t s_insert_seq = seq_suffix[0];
      
      RadixNodeUptr<N> inserted_node(new RadixNode);
      inserted_node->child_nodes[s_insert_branch] = std::move(node->child_nodes[s]);
      inserted_node->child_nodes[s_insert_branch]->branch = std::move(branch_suffix);
      inserted_node->child_nodes[s_insert_seq].reset(new RadixNode);
      inserted_node->child_nodes[s_insert_seq]->branch = std::move(seq_suffix);
      set_insert<T>(inserted_node->child_nodes[s_insert_seq]->terminal_idx, seq_idx);
      node->child_nodes[s] = std::move(inserted_node);
      node->child_nodes[s]->branch = std::move(branch_prefix);
    }
  }
  struct Levenshtein {
    RadixNodeUptr<N> & root;
    ustring_view sequence;
    const int max_distance;
    std::vector<std::pair<uint32_t, int>> output;
    Levenshtein(RadixNodeUptr<N> & root, ustring_view sequence, const int max_distance) : root(root), sequence(sequence), max_distance(max_distance) {}
    std::vector<std::pair<uint32_t, int>> search() {
      auto start_row = int_range(0, sequence.size());
      search(root, start_row);
      return std::move(output); // use std::move, output is not local and RVO doesn't apply
    }
    void search(RadixNodeUptr<N> & node, std::vector<int> & previous_row) {
      if( *std::min_element(previous_row.begin(), previous_row.end()) > max_distance ) {
        return;
      }
      for(auto i : node->terminal_idx) {
        if((previous_row.back() <= max_distance)) {
          output.push_back(std::make_pair(i, previous_row.back()));
        }
      }
      for (size_t i=0; i<N; ++i) {
        if(! node->child_nodes[i]) continue;
        std::vector<int> current_row(sequence.size() + 1);
        ustring_view branch = node->child_nodes[i]->branch;
        update_row(branch[0], previous_row, current_row);
        for(size_t u=1; u<branch.size(); ++u) {
          update_row(branch[u], current_row);
        }
        search(node->child_nodes[i], current_row);
      }
    }
    inline void update_row(const uint8_t branchval, const std::vector<int> & previous_row, std::vector<int> & current_row) {
      current_row[0] = previous_row[0] + 1;
      for(size_t i=1; i<current_row.size(); ++i) {
        int match_cost  = previous_row[i-1] + (sequence[i-1] == branchval ? 0 : 1);
        int insert_cost = current_row[i-1] + 1;
        int delete_cost = previous_row[i] + 1;
        current_row[i] = std::min({match_cost, insert_cost, delete_cost});
      }
    }
    inline void update_row(const uint8_t branchval, std::vector<int> & row) {
      int previous_row_i_minus_1 = row[0];
      row[0] = row[0] + 1;
      for(size_t i=1; i<row.size(); ++i) {
        int match_cost  = previous_row_i_minus_1 + (sequence[i-1] == branchval ? 0 : 1);
        int insert_cost = row[i-1] + 1;
        int delete_cost = row[i] + 1;
        previous_row_i_minus_1 = row[i];
        row[i] = std::min({match_cost, insert_cost, delete_cost});
      }
    }
  };
  
};


namespace DNATree {
const std::string alphabet = "ACGT";
constexpr size_t N = 4;

inline ustring convertString(const std::string & sequence) {
  ustring x(sequence.size(), 0);
  for(size_t j=0; j<sequence.size(); ++j) {
    switch(sequence[j]) {
    case 'A':
      x[j] = 0;
      break;
    case 'C':
      x[j] = 1;
      break;
    case 'G':
      x[j] = 2;
      break;
    case 'T':
      x[j] = 3;
      break;
    default:
      throw std::runtime_error("sequence must have only A, C, G or T");
    }
  }
  return x;
}
inline std::string unconvertString(const ustring & sequence) {
  std::string x(sequence.size(), 0);
  for(size_t j=0; j<sequence.size(); ++j) {
    x[j] = alphabet[sequence[j]];
  }
  return x;
}

std::vector<ustring> convertStrings(const std::vector<std::string> & sequences) {
  std::vector<ustring> output(sequences.size());
  for(size_t i=0; i<sequences.size(); ++i) {
    output[i] = convertString(sequences[i]);
  }
  return output;
}

std::vector<std::string> unconvertStrings(const std::vector<ustring> & sequences) {
  std::vector<std::string> output(sequences.size());
  for(size_t i=0; i<sequences.size(); ++i) {
    output[i] = unconvertString(sequences[i]);
  }
  return output;
}

void print(RadixNodeUptr<N> & node, size_t depth = 0) {
  if(depth == 0) {
    std::cout << "(root)";
  } else {
    for(uint8_t i=0; i<depth; ++i) std::cout << " ";
  }
  std::cout << unconvertString(node->branch) << ":";
  for(auto it = node->terminal_idx.begin(); it != node->terminal_idx.end(); ++it) {
    if(it != node->terminal_idx.begin()) {
      std::cout << ",";
    }
    std::cout << *it;
  }
  std::cout << std::endl;
  for(uint8_t i=0; i<N; ++i) {
    if(node->child_nodes[i]) {
      print(node->child_nodes[i], depth + 1);
    }
  }
}
}

// [[Rcpp::export(rng = false)]]
void og_printTree(const std::vector<std::string> & sequences) {
  RadixNodeUptr<4> root(new RadixNode<4>);
  auto sequences2 = DNATree::convertStrings(sequences);
  for(uint32_t i=0; i<sequences.size(); ++i) {
    RadixNode<4>::insert(root, sequences2[i], i);
  }
  DNATree::print(root);
}


// [[Rcpp::export(rng = false)]]
SEXP og_levenshteinSearch(const std::vector<std::string> & query, const std::vector<std::string> & target, const int max_distance) {
  RadixNodeUptr<4> root(new RadixNode<4>);
  auto query2 = DNATree::convertStrings(query);
  auto target2 = DNATree::convertStrings(target);
  for(uint32_t i=0; i<target2.size(); ++i) {
    RadixNode<4>::insert(root, target2[i], i);
  }
  std::vector<uint32_t> query_idx;
  std::vector<uint32_t> target_idx;
  std::vector<int> distance;
  for(uint32_t i=0; i<query2.size(); ++i) {
    auto res = RadixNode<4>::Levenshtein{root, query2[i], max_distance}.search();
    for(auto & r : res) {
      query_idx.push_back(i);
      target_idx.push_back(r.first);
      distance.push_back(r.second);
    }
  }
  if(query_idx.size() == 0) return(R_NilValue);
  return DataFrame::create(_["query"] = Rcpp::wrap(query_idx), _["target"] = Rcpp::wrap(target_idx), _["distance"] = Rcpp::wrap(distance), _["stringsAsFactors"] = false);
}


// [[Rcpp::export(rng = false)]]
int og_treeBuild(const std::vector<std::string> & target, const int max_distance) {
  RadixNodeUptr<4> root(new RadixNode<4>);
  auto target2 = DNATree::convertStrings(target);
  for(uint32_t i=0; i<target2.size(); ++i) {
    RadixNode<4>::insert(root, target2[i], i);
  }
  return root->size();
}
