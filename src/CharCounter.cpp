#include "seqtrie_types.h"

// [[Rcpp::export(rng = false)]]
CharCounterXPtr CharCounter_create() {
  return CharCounterXPtr(new CharCounter(), true);
}
// [[Rcpp::export(rng = false)]]
void CharCounter_add(CharCounterXPtr xp, CharacterVector sequences) {
  auto & counts = *xp;
  size_t nseqs = Rf_xlength(sequences);
  const SEXP * seqs = STRING_PTR_RO(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    auto s = charsxp_to_cspan(seqs[i]);
    for(auto c : s) {
      counts[c]++; // value initialized to zero if not present
    }
  }
}
// [[Rcpp::export(rng = false)]]
void CharCounter_subtract(CharCounterXPtr xp, CharacterVector sequences) {
  auto & counts = *xp;
  size_t nseqs = Rf_xlength(sequences);
  const SEXP * seqs = STRING_PTR_RO(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    auto s = charsxp_to_cspan(seqs[i]);
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
}
// [[Rcpp::export(rng = false)]]
CharacterVector CharCounter_keys(CharCounterXPtr xp) {
  auto & counts = *xp;
  CharacterVector out(counts.size());
  size_t i=0;
  for(auto & kv : counts) {
    out[i] = Rf_mkCharLen(&kv.first, 1);
    i++;
  }
  return out;
}
// [[Rcpp::export(rng = false)]]
CharacterVector get_charset(CharacterVector sequences) {
  size_t nseqs = Rf_xlength(sequences);
  std::set<char> charset;
  const SEXP * seqs = STRING_PTR_RO(sequences);
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
