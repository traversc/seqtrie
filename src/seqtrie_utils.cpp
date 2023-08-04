#include <Rcpp.h>
#include "seqtrie_types.h"

using namespace Rcpp;

// Convert a string to a SEXP
// Be careful about R protection / GC
SEXP to_charsxp(const trqwe::small_array<char> & x) {
    return Rf_mkCharLen(x.data(), x.size());
}

// Define a span of const char from a SEXP
cspan charsxp_to_cspan(SEXP x) {
    return cspan(CHAR(x), Rf_xlength(x));
}

// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map
pairchar_map convert_cost_matrix(IntegerMatrix cost_matrix) {
  pairchar_map cost_map;
  std::vector<char> map_elements;
  {
    List dimnames = cost_matrix.attr("dimnames");
    CharacterVector rownames = dimnames[0];
    map_elements.resize(rownames.size());
    for(size_t i=0; i<map_elements.size(); ++i) {
      if(rownames[i] == "gap") {
        map_elements[i] = '\0';
      } else {
        Rcpp::String s = rownames[i];
        map_elements[i] = s.get_cstring()[0];
      }
    }
  }
  size_t N = map_elements.size();
  for(int i=0; i<N; ++i) {
    for(int j=0; j<N; ++j) {
      cost_map[pairchar(map_elements[i], map_elements[j])] = cost_matrix(i,j);
    }
  }
  return cost_map;
}
