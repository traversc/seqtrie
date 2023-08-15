#include <set>
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
  for(size_t i=0; i<N; ++i) {
    for(size_t j=0; j<N; ++j) {
      cost_map[pairchar(map_elements[i], map_elements[j])] = cost_matrix(i,j);
    }
  }
  return cost_map;
}

void CharCounter::add(cspan s) {
  for(auto c : s) {
    counts[c]++; // value initialized to zero if not present
  }
}
void CharCounter::subtract(cspan s) {
  for(auto c : s) {
    size_t & cc = counts[c];
    if(cc > 0) {
      cc--;
    }
    if(cc == 0) { // if the count becomes zero, remove the key
      counts.erase(c);
    }
  }
}
CharacterVector CharCounter::keys() {
  CharacterVector out(counts.size());
  size_t i=0;
  for(auto & kv : counts) {
    out[i] = Rf_mkCharLen(&kv.first, 1);
    i++;
  }
  return out;
}
// [[Rcpp::export(rng = false)]]
CharCounterXPtr CharCounter_create() {
  return CharCounterXPtr(new CharCounter(), true);
}
// [[Rcpp::export(rng = false)]]
void CharCounter_add(CharCounterXPtr xp, CharacterVector sequences) {
  auto & x = *xp;
  size_t nseqs = Rf_xlength(sequences);
  SEXP * seqs = STRING_PTR(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    x.add(charsxp_to_cspan(seqs[i]));
  }
}
// [[Rcpp::export(rng = false)]]
void CharCounter_subtract(CharCounterXPtr xp, CharacterVector sequences) {
  auto & x = *xp;
  size_t nseqs = Rf_xlength(sequences);
  SEXP * seqs = STRING_PTR(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    x.subtract(charsxp_to_cspan(seqs[i]));
  }
}
// [[Rcpp::export(rng = false)]]
CharacterVector CharCounter_keys(CharCounterXPtr xp) {
  return xp->keys();
}
// [[Rcpp::export(rng = false)]]
CharacterVector get_charset(CharacterVector sequences) {
  size_t nseqs = Rf_xlength(sequences);
  std::set<char> charset;
  SEXP * seqs = STRING_PTR(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    cspan x = charsxp_to_cspan(seqs[i]);
    charset.insert(x.begin(), x.end());
  }
  CharacterVector out(charset.size());
  size_t i = 0;
  for(auto c : charset) {
    SET_STRING_ELT(out, i, Rf_mkCharLen(&c, 1));
    i++;
  }
  return out;
}