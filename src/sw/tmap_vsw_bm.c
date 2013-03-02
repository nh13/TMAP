/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <config.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../util/tmap_sam_convert.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sort.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_rand.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../sw/tmap_sw.h"
#include "../sw/tmap_fsw.h"
#include "../sw/tmap_vsw.h"
#include "../sw/ksw.h"
#include "../map/util/tmap_map_opt.h"
#include "../map/util/tmap_map_util.h"

#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
void
tmap_vsw_bm_core(int32_t seq_len, int32_t tlen, int32_t n_iter,
                 int32_t n_sub_iter, int32_t vsw_type, int32_t bw)
{
  int32_t i, j, k;
  tmap_vsw_t *vsw = NULL;
  tmap_vsw_opt_t *vsw_opt = NULL;
  int32_t softclip_start, softclip_end;
  tmap_sw_param_t ap;
  int32_t matrix32[25];
  int8_t matrix8[25];

  tmap_map_opt_t *opt = tmap_map_opt_init(TMAP_MAP_ALGO_NONE);

  uint8_t *seq, *target;
  tmap_rand_t *rand = tmap_rand_init(0);

  seq = tmap_malloc(sizeof(uint8_t) * seq_len, "seq");
  target = tmap_malloc(sizeof(uint8_t) * tlen, "target");

  // random sequence
  for(i=0;i<seq_len;i++) {
      seq[i] = (uint8_t)(4*tmap_rand_get(rand));
  }

  softclip_start = 1;
  softclip_end = 1;

  // initialize opt
  if(0 <= vsw_type) { 
      vsw_opt = tmap_vsw_opt_init(opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape, opt->score_thr);
      vsw = tmap_vsw_init(seq, seq_len, softclip_start, softclip_end, vsw_type, vsw_opt);
  }
  else {
      ap.matrix = matrix32;
      __map_util_gen_ap(ap, opt);
      for(i=0;i<5;i++) {
          for(j=0;j<5;j++) {
              if(i==j) matrix8[i + (i*j)] = opt->score_match;
              else matrix8[i + (i*j)] = opt->pen_mm;
          }
      }
  }


  int32_t front = (tlen - seq_len) / 2;
  int32_t end = tlen - seq_len - front;
  while(i<n_iter) {
      tmap_map_sam_t tmp_sam;
      int32_t overflow;
      for(j=k=0;j<front;j++,k++) {
          target[k] = (uint8_t)(4*tmap_rand_get(rand));
      }
      for(j=0;j<seq_len;j++,k++) {
          target[k] = seq[j];
      }
      for(j=0;j<end;j++,k++) {
          target[k] = (uint8_t)(4*tmap_rand_get(rand));
      }
      for(j=0;j<n_sub_iter&&i<n_iter;j++,i++) {
          if(0 <= vsw_type) { 
              // initialize the bounds
              tmp_sam.result.query_start = tmp_sam.result.query_end = 0;
              tmp_sam.result.target_start = tmp_sam.result.target_end = 0;
              // run the vsw
              tmap_vsw_process_fwd(vsw, seq, seq_len, target, tlen,
                            &tmp_sam.result, &overflow, opt->score_thr, 0);
          }
          else if(-1 == vsw_type) {
              tmap_sw_clipping_core(seq, seq_len, target, tlen,
                                    &ap, softclip_start, softclip_end,
                                    NULL, NULL, 0);
          }
          else if(-2 == vsw_type) {
              ksw_global(seq_len, seq, tlen, target, 4, matrix8, opt->pen_gapo, opt->pen_gape, bw, NULL, NULL); 
          }
      }
  }

  // free memory
  free(target);
  free(seq);
  if(0 <= vsw_type) {
      tmap_vsw_opt_destroy(vsw_opt);
      tmap_vsw_destroy(vsw);
  }
  tmap_map_opt_destroy(opt);
  tmap_rand_destroy(rand);
}

static int
usage(int32_t seq_len, int32_t tlen, int32_t n_iter, 
      int32_t n_sub_iter, int32_t vsw_type)
{
  static char *vsw_type_str[] = {
      "NB: currently only #1, #4, and #6 have been tested",
      "-2 - ksw_global",
      "-1 - tmap_sw_clipping_core",
      "1 - lh3/ksw.c/nh13",
      "2 - simple VSW",
      "3 - SHRiMP2 VSW [not working]",
      "4 - Psyho (Top Coder #1)",
      "5 - ACRush (Top Coder #2)",
      "6 - folsena (Top Coder #3)",
      "7 - logicmachine (Top Coder #4)",
      "8 - venco (Top Coder #5) [not working]",
      "9 - Bladze (Top Coder #6)",
      "10 - ngthuydiem (Top Coder #7) [Farrar cut-and-paste]",
      NULL};
  char **ptr = NULL;
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s vswbm [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -q INT      the query length [%d]\n", seq_len);
  tmap_file_fprintf(tmap_file_stderr, "         -t INT      the target length [%d] (must be at least as long as the query)\n", tlen);
  tmap_file_fprintf(tmap_file_stderr, "         -n INT      the number of iterations [%d]\n", n_iter);
  tmap_file_fprintf(tmap_file_stderr, "         -N INT      the number of re-evaluations of the same query/target combination [%d]\n", n_sub_iter);
  tmap_file_fprintf(tmap_file_stderr, "         -H INT      smith waterman algorithm [%d]\n", vsw_type);
  ptr = vsw_type_str;
  while(NULL != (*ptr)) {
      tmap_file_fprintf(tmap_file_stderr, "                       %s\n", (*ptr));
      ptr++;
  }
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int
tmap_vswbm_main(int argc, char *argv[])
{

  int32_t seq_len = 150;
  int32_t tlen = 256; 
  int32_t n_iter = 1000;
  int32_t n_sub_iter = 1;
  int32_t vsw_type = 1;
  int32_t bw = 50;
  int c;

  while((c = getopt(argc, argv, "q:t:n:N:H:w:h")) >= 0) {
      switch(c) {
        case 'q':
          seq_len = atoi(optarg); break;
        case 't':
          tlen = atoi(optarg); break;
        case 'n':
          n_iter = atoi(optarg); break;
        case 'N':
          n_sub_iter = atoi(optarg); break;
        case 'H':
          vsw_type = atoi(optarg); break;
        case 'w':
          bw = atoi(optarg); break;
        case 'h':
        default:
          return usage(seq_len, tlen, n_iter, n_sub_iter, vsw_type);
      }
  }
     
  if(argc != optind || seq_len > tlen) {
      return usage(seq_len, tlen, n_iter, n_sub_iter, vsw_type);
  }

  tmap_progress_set_verbosity(1);
  tmap_progress_print2("starting benchmark");

  tmap_vsw_bm_core(seq_len, tlen, n_iter, n_sub_iter, vsw_type, bw);
  
  tmap_progress_print2("ending benchmark");

  return 0;
}
#endif
