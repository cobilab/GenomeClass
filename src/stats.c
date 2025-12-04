#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "defs.h"
#include "io.h"
#include "misc.h"
#include "stats.h"
#include "auxgenomeclass.h"

#define MAX_SYMBOLS   65536u
#define BIGRAM_SIZE_8 (256u * 256u)

/* --------------------------------------------------------------------
 * Hash table for 16-bit bigrams
 * ------------------------------------------------------------------*/

typedef struct
  {
  uint32_t key;   // (symbol1 << 16) | symbol2
  uint64_t count;
  int      used;
  }
BigramEntry16;

typedef struct
  {
  BigramEntry16 *table;
  size_t         capacity;
  size_t         size;
  }
BigramHash16;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint32_t hash32(uint32_t x)
  {
  x ^= x >> 16;
  x *= 0x7feb352d;
  x ^= x >> 15;
  x *= 0x846ca68b;
  x ^= x >> 16;
  return x;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static size_t next_power_of_two(size_t x)
  {
  size_t p = 1;
  while(p < x)
    p <<= 1;
  return p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void bigram_hash_init(BigramHash16 *h, size_t initial_capacity)
  {
  size_t cap = next_power_of_two(initial_capacity);
  if(cap < 1024)
    cap = 1024;

  h->capacity = cap;
  h->size     = 0;
  h->table    = (BigramEntry16 *)calloc(h->capacity, sizeof(BigramEntry16));
  if(!h->table)
    {
    fprintf(stderr, "Memory allocation failed in bigram_hash_init\n");
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void bigram_hash_free(BigramHash16 *h)
  {
  free(h->table);
  h->table    = NULL;
  h->capacity = 0;
  h->size     = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void bigram_hash_insert_into_table(
  BigramEntry16 *table,
  size_t capacity,
  uint32_t key,
  uint64_t count
  )
  {
  uint32_t h = hash32(key);
  size_t idx = (size_t)(h & (uint32_t)(capacity - 1));

  for(;;)
    {
    if(!table[idx].used)
      {
      table[idx].used  = 1;
      table[idx].key   = key;
      table[idx].count = count;
      return;
      }
    if(table[idx].key == key)
      {
      table[idx].count += count;
      return;
      }
    idx = (idx + 1) & (capacity - 1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void bigram_hash_rehash(BigramHash16 *h, size_t new_capacity)
  {
  BigramEntry16 *new_table;
  size_t cap = next_power_of_two(new_capacity);

  new_table = (BigramEntry16 *)calloc(cap, sizeof(BigramEntry16));
  if(!new_table)
    {
    fprintf(stderr, "Memory allocation failed in bigram_hash_rehash\n");
    exit(1);
    }

  for(size_t i = 0 ; i < h->capacity ; ++i)
    {
    if(h->table[i].used)
      bigram_hash_insert_into_table(new_table, cap,
                                    h->table[i].key,
                                    h->table[i].count);
    }

  free(h->table);
  h->table    = new_table;
  h->capacity = cap;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void bigram_hash_inc(BigramHash16 *h, uint32_t key)
  {
  /* Resize when load factor > 0.7 */
  if(h->size * 10 >= h->capacity * 7)
    bigram_hash_rehash(h, h->capacity * 2);

  uint32_t hv = hash32(key);
  size_t idx = (size_t)(hv & (uint32_t)(h->capacity - 1));

  for(;;)
    {
    if(!h->table[idx].used)
      {
      h->table[idx].used  = 1;
      h->table[idx].key   = key;
      h->table[idx].count = 1;
      h->size++;
      return;
      }
    if(h->table[idx].key == key)
      {
      h->table[idx].count++;
      return;
      }
    idx = (idx + 1) & (h->capacity - 1);
    }
  }

/* --------------------------------------------------------------------
 * Single-pass readers for 8-bit and 16-bit symbols
 * ------------------------------------------------------------------*/

static int read_file_and_collect_8(
  const char   *filename,
  uint64_t      counts[MAX_SYMBOLS],
  uint64_t      bigram_counts[BIGRAM_SIZE_8],
  uint64_t     *N_out,
  uint64_t     *run_count_out,
  uint64_t     *max_run_len_out,
  uint64_t     *num_changes_out,
  uint64_t     *equal_pairs_out,
  uint64_t     *bigram_total_out
  )
  {
  FILE *fp = Fopen(filename, "rb");
  if(!fp)
    {
    perror("Error opening file");
    return 1;
    }

  uint64_t N            = 0;
  uint64_t run_count    = 0;
  uint64_t max_run_len  = 0;
  uint64_t num_changes  = 0;
  uint64_t equal_pairs  = 0;
  uint64_t bigram_total = 0;

  int have_prev         = 0;
  unsigned char prev    = 0;
  uint64_t current_run_len = 0;

  unsigned char buffer[BUFFER_SIZE];
  size_t nread;

  while((nread = fread(buffer, 1, BUFFER_SIZE, fp)) > 0)
    {
    for(size_t i = 0 ; i < nread ; ++i)
      {
      unsigned char c = buffer[i];

      /* Count current symbol */
      counts[c]++;
      N++;

      if(!have_prev)
        {
        /* First symbol */
        prev = c;
        have_prev = 1;
        current_run_len = 1;
        run_count = 1;
        }
      else
        {
        /* Runs and changes */
        if(c == prev)
          {
          current_run_len++;
          equal_pairs++;
          }
        else
          {
          num_changes++;
          if(current_run_len > max_run_len)
            max_run_len = current_run_len;
          current_run_len = 1;
          run_count++;
          }

        /* Bigrams: (prev, c) */
        unsigned int idx = ((unsigned int)prev << 8) | c;
        bigram_counts[idx]++;
        bigram_total++;

        prev = c;
        }
      }
    }

  fclose(fp);

  if(have_prev && current_run_len > max_run_len)
    max_run_len = current_run_len;

  *N_out            = N;
  *run_count_out    = run_count;
  *max_run_len_out  = max_run_len;
  *num_changes_out  = num_changes;
  *equal_pairs_out  = equal_pairs;
  *bigram_total_out = bigram_total;

  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static int read_file_and_collect_16(
  const char   *filename,
  uint64_t      counts[MAX_SYMBOLS],
  BigramHash16 *bh,
  uint64_t     *N_out,
  uint64_t     *run_count_out,
  uint64_t     *max_run_len_out,
  uint64_t     *num_changes_out,
  uint64_t     *equal_pairs_out,
  uint64_t     *bigram_total_out
  )
  {
  FILE *fp = Fopen(filename, "rb");
  if(!fp)
    {
    perror("Error opening file");
    return 1;
    }

  uint64_t N            = 0;
  uint64_t run_count    = 0;
  uint64_t max_run_len  = 0;
  uint64_t num_changes  = 0;
  uint64_t equal_pairs  = 0;
  uint64_t bigram_total = 0;

  int have_prev         = 0;
  uint16_t prev         = 0;
  uint64_t current_run_len = 0;

  int have_half         = 0;
  unsigned char hi_byte = 0;

  unsigned char buffer[BUFFER_SIZE];
  size_t nread;

  while((nread = fread(buffer, 1, BUFFER_SIZE, fp)) > 0)
    {
    for(size_t i = 0 ; i < nread ; ++i)
      {
      unsigned char b = buffer[i];

      if(!have_half)
        {
        hi_byte = b;
        have_half = 1;
        }
      else
        {
        uint16_t sym = ((uint16_t)hi_byte << 8) | b;
        have_half = 0;

        counts[sym]++;
        N++;

        if(!have_prev)
          {
          prev = sym;
          have_prev = 1;
          current_run_len = 1;
          run_count = 1;
          }
        else
          {
          if(sym == prev)
            {
            current_run_len++;
            equal_pairs++;
            }
          else
            {
            num_changes++;
            if(current_run_len > max_run_len)
              max_run_len = current_run_len;
            current_run_len = 1;
            run_count++;
            }

          /* Bigram key: (prev, sym) */
          uint32_t key = ((uint32_t)prev << 16) | (uint32_t)sym;
          bigram_hash_inc(bh, key);
          bigram_total++;

          prev = sym;
          }
        }
      }
    }

  fclose(fp);

  if(have_prev && current_run_len > max_run_len)
    max_run_len = current_run_len;

  *N_out            = N;
  *run_count_out    = run_count;
  *max_run_len_out  = max_run_len;
  *num_changes_out  = num_changes;
  *equal_pairs_out  = equal_pairs;
  *bigram_total_out = bigram_total;

  return 0;
  }

/* --------------------------------------------------------------------
 * Basic stats, run stats, bigram stats
 * ------------------------------------------------------------------*/

static void compute_basic_stats(
  const uint64_t counts[MAX_SYMBOLS],
  uint32_t alphabet_size,
  uint64_t N,
  BasicStats *out
  )
  {
  uint64_t V = 0;
  uint64_t max_freq = 0;
  uint64_t min_freq_nonzero = 0;

  for(uint32_t i = 0 ; i < alphabet_size ; ++i)
    {
    if(counts[i] > 0)
      {
      V++;
      if(counts[i] > max_freq)
        max_freq = counts[i];
      if(min_freq_nonzero == 0 || counts[i] < min_freq_nonzero)
        min_freq_nonzero = counts[i];
      }
    }

  out->N                = N;
  out->V                = V;
  out->max_freq         = max_freq;
  out->min_freq_nonzero = min_freq_nonzero;
  out->diversity_rel    = 0.0;
  out->freq_rel_max     = 0.0;
  out->diff_max_min_rel = 0.0;

  if(N > 0)
    {
    out->diversity_rel = (double)V / (double)N;
    out->freq_rel_max  = (double)max_freq / (double)N;
    if(min_freq_nonzero > 0)
      out->diff_max_min_rel =
        (double)(max_freq - min_freq_nonzero) / (double)N;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void compute_run_stats(
  uint64_t N,
  uint64_t run_count,
  uint64_t max_run_len,
  uint64_t num_changes,
  RunStats *out
  )
  {
  out->max_run_len = max_run_len;
  out->run_count   = run_count;
  out->num_changes = num_changes;
  out->avg_run_len = 0.0;
  out->change_rate = 0.0;

  if(run_count > 0)
    out->avg_run_len = (double)N / (double)run_count;

  if(N > 1)
    out->change_rate = (double)num_changes / (double)(N - 1);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void compute_bigram_stats_8(
  const uint64_t bigram_counts[BIGRAM_SIZE_8],
  uint64_t N,
  uint64_t equal_pairs,
  uint64_t bigram_total,
  BigramStats *out
  )
  {
  out->B               = 0;
  out->equal_pairs     = equal_pairs;
  out->bigram_total    = bigram_total;
  out->equal_pair_ratio = 0.0;
  out->bigram_div_rel   = 0.0;

  for(int k = 0 ; k < 3 ; ++k)
    {
    out->top_count[k] = 0;
    out->top_key[k]   = 0;
    }

  if(N > 1)
    out->equal_pair_ratio = (double)equal_pairs / (double)(N - 1);

  for(int idx = 0 ; idx < (int)BIGRAM_SIZE_8 ; ++idx)
    {
    uint64_t c = bigram_counts[idx];
    if(c > 0)
      {
      out->B++;

      uint32_t s1  = (uint32_t)(idx >> 8);
      uint32_t s2  = (uint32_t)(idx & 0xFFu);
      uint32_t key = (s1 << 16) | s2;

      if(c > out->top_count[0])
        {
        out->top_count[2] = out->top_count[1];
        out->top_key[2]   = out->top_key[1];
        out->top_count[1] = out->top_count[0];
        out->top_key[1]   = out->top_key[0];
        out->top_count[0] = c;
        out->top_key[0]   = key;
        }
      else if(c > out->top_count[1])
        {
        out->top_count[2] = out->top_count[1];
        out->top_key[2]   = out->top_key[1];
        out->top_count[1] = c;
        out->top_key[1]   = key;
        }
      else if(c > out->top_count[2])
        {
        out->top_count[2] = c;
        out->top_key[2]   = key;
        }
      }
    }

  if(bigram_total > 0)
    out->bigram_div_rel = (double)out->B / (double)bigram_total;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void compute_bigram_stats_16(
  const BigramHash16 *bh,
  uint64_t N,
  uint64_t equal_pairs,
  uint64_t bigram_total,
  BigramStats *out
  )
  {
  out->B               = bh->size;
  out->equal_pairs     = equal_pairs;
  out->bigram_total    = bigram_total;
  out->equal_pair_ratio = 0.0;
  out->bigram_div_rel   = 0.0;

  for(int k = 0 ; k < 3 ; ++k)
    {
    out->top_count[k] = 0;
    out->top_key[k]   = 0;
    }

  if(N > 1)
    out->equal_pair_ratio = (double)equal_pairs / (double)(N - 1);

  for(size_t i = 0 ; i < bh->capacity ; ++i)
    {
    if(!bh->table[i].used)
      continue;

    uint64_t c   = bh->table[i].count;
    uint32_t key = bh->table[i].key;

    if(c > out->top_count[0])
      {
      out->top_count[2] = out->top_count[1];
      out->top_key[2]   = out->top_key[1];
      out->top_count[1] = out->top_count[0];
      out->top_key[1]   = out->top_key[0];
      out->top_count[0] = c;
      out->top_key[0]   = key;
      }
    else if(c > out->top_count[1])
      {
      out->top_count[2] = out->top_count[1];
      out->top_key[2]   = out->top_key[1];
      out->top_count[1] = c;
      out->top_key[1]   = key;
      }
    else if(c > out->top_count[2])
      {
      out->top_count[2] = c;
      out->top_key[2]   = key;
      }
    }

  if(bigram_total > 0)
    out->bigram_div_rel = (double)out->B / (double)bigram_total;
  }


/* --------------------------------------------------------------------
 * Usage for the stats subcommand
 * ------------------------------------------------------------------*/

static void print_stats_usage(const char *progname)
  {
  printf("Usage: %s stats [-b 8|16] <filename>\n", progname);
  printf("  -b 8|16     Symbol width (8-bit or 16-bit)\n");
  }

/* --------------------------------------------------------------------
 * Entry point used from ox.c
 * ------------------------------------------------------------------*/

char * stats_mode(char *filename, int bits)
  {
  int opt;

  static uint64_t counts[MAX_SYMBOLS];
  static uint64_t bigram_counts8[BIGRAM_SIZE_8];

  memset(counts,        0, sizeof(counts));
  memset(bigram_counts8, 0, sizeof(bigram_counts8));

  uint64_t N, run_count, max_run_len, num_changes;
  uint64_t equal_pairs, bigram_total;

  BasicStats  basic;
  RunStats    runs;
  BigramStats bigrams;

  if(bits == 8)
    {
    if(read_file_and_collect_8(filename,
                               counts,
                               bigram_counts8,
                               &N,
                               &run_count,
                               &max_run_len,
                               &num_changes,
                               &equal_pairs,
                               &bigram_total) != 0)
      exit(1);

    compute_basic_stats(counts, 256u, N, &basic);
    compute_run_stats(N, run_count, max_run_len, num_changes, &runs);
    compute_bigram_stats_8(bigram_counts8,
                           N,
                           equal_pairs,
                           bigram_total,
                           &bigrams);
    }
  else  /* bits == 16 */
    {
    BigramHash16 bh;
    bigram_hash_init(&bh, 1024);

    if(read_file_and_collect_16(filename,
                                counts,
                                &bh,
                                &N,
                                &run_count,
                                &max_run_len,
                                &num_changes,
                                &equal_pairs,
                                &bigram_total) != 0)
      {
      bigram_hash_free(&bh);
      exit(1);
      }

    compute_basic_stats(counts, 65536u, N, &basic);
    compute_run_stats(N, run_count, max_run_len, num_changes, &runs);
    compute_bigram_stats_16(&bh,
                            N,
                            equal_pairs,
                            bigram_total,
                            &bigrams);

    bigram_hash_free(&bh);
    }

  //printf("File for stats: %s\n", filename);
  //printf("Considering Symbol width %d bits\n\n", bits);

  char *results = strdup("");

  //Distinct_symbols\tRelative_diversity\tMost_common_symbol_freq_rel\tMax_min_freq_difference_rel\t

  results = concatenate_strings(results, int_to_string((unsigned long long)basic.V), 1);
  results = concatenate_strings(results, float_to_string(basic.diversity_rel), 1);
  results = concatenate_strings(results, float_to_string(basic.freq_rel_max), 1);


  if(basic.min_freq_nonzero > 0)
    {
    results = concatenate_strings(results, float_to_string(basic.diff_max_min_rel), 1);
    }
  else
    {
    results = concatenate_strings(results, float_to_string(-1.0), 1);
    }


  //Longest_run_length\tAverage_run_length\tChange_rate\t
  results = concatenate_strings(results, int_to_string((unsigned long long)runs.max_run_len), 1);
  results = concatenate_strings(results, float_to_string(runs.avg_run_len), 1);
  results = concatenate_strings(results, float_to_string(runs.change_rate), 1);


  //Equal_adjacent_pairs_ratio\tDistinct_bigrams\tBigram_diversity\tMost_frequent_bigram_1\tFreq_bigram_1\tMost_frequent_bigram_2\tFreq_bigram_2\tMost_frequent_bigram_3\tFreq_bigram_3

  results = concatenate_strings(results, float_to_string(bigrams.equal_pair_ratio), 1);
  results = concatenate_strings(results, int_to_string((unsigned long long)bigrams.B), 1);
  results = concatenate_strings(results, float_to_string(bigrams.bigram_div_rel), 1);

  for(int k = 0 ; k < 3 ; ++k)
    {
    if(bigrams.top_count[k] == 0)
      continue;

    uint32_t key = bigrams.top_key[k];
    uint32_t s1  = (key >> 16) & 0xFFFFu;
    uint32_t s2  = key & 0xFFFFu;
    double rel   = 0.0;
    if(bigrams.bigram_total > 0)
      rel = (double)bigrams.top_count[k] /
            (double)bigrams.bigram_total;

    char buffer[100];
    sprintf(buffer, "(%u, %u)", s1, s2);

    results = concatenate_strings(results, buffer, 1);
    results = concatenate_strings(results, int_to_string((unsigned long long)bigrams.top_count[k]), 1);
    }
  

  return results;
  }


