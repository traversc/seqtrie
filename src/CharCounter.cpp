#include "seqtrie_types.h"

// [[Rcpp::export(rng = false)]]
CharCounterXPtr CharCounter_create() {
  auto ptr = std::make_unique<CharCounter>();
  return CharCounterXPtr(ptr.release(), true);
}
// [[Rcpp::export(rng = false)]]
void CharCounter_add(CharCounterXPtr xp, CharacterVector sequences) {
  auto & counts = *xp;
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  const SEXP * seqs = STRING_PTR_RO(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    if(seqs[i] == NA_STRING) continue;
    auto s = charsxp_to_cspan(seqs[i]);
    for(char c : s) {
      counts.add(c);
    }
  }
}
// [[Rcpp::export(rng = false)]]
void CharCounter_subtract(CharCounterXPtr xp, CharacterVector sequences) {
  auto & counts = *xp;
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  const SEXP * seqs = STRING_PTR_RO(sequences);
  for(size_t i=0; i<nseqs; ++i) {
    if(seqs[i] == NA_STRING) continue;
    auto s = charsxp_to_cspan(seqs[i]);
    for(char c : s) {
      counts.subtract(c);
    }
  }
}
// [[Rcpp::export(rng = false)]]
CharacterVector CharCounter_keys(CharCounterXPtr xp) {
  auto & counts = *xp;
  CharacterVector out(counts.size());
  size_t i = 0;
  for(size_t byte = 0; byte < CharCounter::alphabet_size; ++byte) {
    if(counts.counts[byte] == 0) continue;
    const char c = static_cast<char>(static_cast<unsigned char>(byte));
    SET_STRING_ELT(out, checked_r_xlen(i), Rf_mkCharLen(&c, 1));
    ++i;
  }
  return out;
}
// [[Rcpp::export(rng = false)]]
CharacterVector get_charset(CharacterVector sequences) {
  CharCounter counts;
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  const SEXP * seqs = STRING_PTR_RO(sequences);
  for(size_t i = 0; i < nseqs; ++i) {
    if(seqs[i] == NA_STRING) continue;
    cspan x = charsxp_to_cspan(seqs[i]);
    for(char c : x) {
      counts.add(c);
    }
  }
  CharacterVector out(counts.size());
  size_t i = 0;
  for(size_t byte = 0; byte < CharCounter::alphabet_size; ++byte) {
    if(counts.counts[byte] == 0) continue;
    const char c = static_cast<char>(static_cast<unsigned char>(byte));
    SET_STRING_ELT(out, checked_r_xlen(i), Rf_mkCharLen(&c, 1));
    ++i;
  }
  return out;
}
