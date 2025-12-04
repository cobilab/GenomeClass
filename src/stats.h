#ifndef STATS_H
#define STATS_H

#include <stdint.h>

/* --------------------------------------------------------------------
 * Structures for storing statistics:
 * (1) Basic counts and ratios
 * (2) Repetition and alternation (runs)
 * (3) Local structure (bigrams)
 * ------------------------------------------------------------------*/

typedef struct
  {
  uint64_t N;                 // total length (number of symbols)
  uint64_t V;                 // number of distinct symbols
  uint64_t max_freq;          // absolute frequency of the most common symbol
  uint64_t min_freq_nonzero;  // minimum non-zero frequency among used symbols
  double   diversity_rel;     // V / N
  double   freq_rel_max;      // max_freq / N
  double   diff_max_min_rel;  // (max_freq - min_freq_nonzero) / N
  }
BasicStats;

typedef struct
  {
  uint64_t max_run_len;   // longest run of equal symbols
  uint64_t run_count;     // number of runs
  uint64_t num_changes;   // number of changes s[i] != s[i+1]
  double   avg_run_len;   // N / run_count
  double   change_rate;   // num_changes / (N-1)
  }
RunStats;

typedef struct
  {
  uint64_t B;                 // number of distinct bigrams
  uint64_t equal_pairs;       // number of pairs with s[i] == s[i+1]
  uint64_t bigram_total;      // total number of bigrams (≈ N-1)
  double   equal_pair_ratio;  // equal_pairs / (N-1)
  double   bigram_div_rel;    // B / bigram_total
  uint64_t top_count[3];      // counts of the top-3 bigrams
  uint32_t top_key[3];        // (symbol1<<16 | symbol2) for the top-3 bigrams
  }
BigramStats;

/* Entry point used by `ox stats ...` */
char* stats_mode(char* filename, int bits);

#endif
