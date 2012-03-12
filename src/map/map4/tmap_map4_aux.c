#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_vec.h"
#include "../../util/tmap_rand.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../index/tmap_bwt_smem.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_util.h"
#include "tmap_map4.h"
#include "tmap_map4_aux.h"

tmap_map4_aux_smem_iter_t *
tmap_map4_aux_smem_iter_init()
{
  tmap_map4_aux_smem_iter_t *iter;
  iter = tmap_calloc(1, sizeof(tmap_map4_aux_smem_iter_t), "iter");
  iter->tmpvec[0] = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->tmpvec[0]");
  iter->tmpvec[1] = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->tmpvec[1]");
  iter->matches   = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->matches");
  return iter;
}

void 
tmap_map4_aux_smem_iter_destroy(tmap_map4_aux_smem_iter_t *iter)
{
  free(iter->tmpvec[0]->a);
  free(iter->tmpvec[1]->a);
  free(iter->tmpvec[0]);
  free(iter->tmpvec[1]);
  free(iter->matches->a);
  free(iter->matches);
  free(iter);
}

static void 
tmap_map4_aux_smem_iter_set_query(tmap_map4_aux_smem_iter_t *iter, int len, const uint8_t *query)
{
  iter->query = query;
  iter->start = 0;
  iter->len = len;
}

static int
tmap_map4_aux_smem_iter_next(tmap_map4_aux_smem_iter_t *iter, tmap_bwt_t *bwt)
{
  iter->tmpvec[0]->n = iter->tmpvec[1]->n = iter->matches->n = 0;
  if (iter->start >= iter->len || iter->start < 0) return -1;
  while (iter->start < iter->len && iter->query[iter->start] > 3) ++iter->start; // skip ambiguous bases
  if (iter->start == iter->len) return -1;
  iter->start = tmap_bwt_smem1(bwt, iter->len, iter->query, iter->start, iter->matches, iter->tmpvec);
  return iter->start;
}

static void
tmap_bwt_smem_intv_copy(tmap_bwt_smem_intv_t *dest, tmap_bwt_smem_intv_t *src)
{
  dest->x[0] = src->x[0];
  dest->x[1] = src->x[1];
  dest->size = src->size;
  dest->info = src->info;
}

static tmap_bwt_smem_intv_vec_t *
tmap_bwt_smem_intv_vec_init()
{
  tmap_bwt_smem_intv_vec_t *matches;
  matches = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "matches");
  matches->a = NULL;
  matches->n = matches->m = 0;
  return matches;
}

static void
tmap_bwt_smem_intv_vec_push(tmap_bwt_smem_intv_vec_t *matches, tmap_bwt_smem_intv_t *p)
{
  if (matches->n == matches->m) {
      matches->m = (0 < matches->m) ? (matches->m << 1) : 2;
      matches->a = tmap_realloc(matches->a, sizeof(tmap_bwt_smem_intv_t) * matches->m, "matches->a");
  }
  tmap_bwt_smem_intv_copy(&matches->a[matches->n++], p);
}

static void
tmap_bwt_smem_intv_vec_destroy(tmap_bwt_smem_intv_vec_t *matches)
{
  free(matches->a);
  free(matches);
}

tmap_map_sams_t *
tmap_map4_aux_core(tmap_seq_t *seq,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_bwt_match_hash_t *hash,
                   tmap_map4_aux_smem_iter_t *iter,
                   tmap_rand_t *rand,
                   tmap_map_opt_t *opt)
{
  int32_t i, n;
  int32_t start, by;
  int32_t min_seed_length, max_seed_length;
  tmap_bwt_int_t k;
  tmap_map_sams_t *sams;
  tmap_bwt_smem_intv_vec_t *matches;
  int32_t total = 0;
  int32_t max_repr;
  
  uint8_t *query;
  int32_t query_len;

  sams = tmap_map_sams_init(NULL);

  query = (uint8_t*)tmap_seq_get_bases(seq)->s;
  query_len = tmap_seq_get_bases_length(seq);

  if(query_len <= 0) return sams;
  
  matches = tmap_bwt_smem_intv_vec_init();

  min_seed_length = opt->min_seed_length;
  max_seed_length = opt->max_seed_length;

  if(query_len < min_seed_length) min_seed_length = query_len; 
  if (min_seed_length < 0) tmap_bug(); // this should not happen, or fix it upstream
  if(-1 == max_seed_length || query_len < max_seed_length) max_seed_length = query_len; 

  start = 0;
  by = (opt->seed_step < 0) ? query_len : opt->seed_step;

  max_repr = opt->max_repr;
  max_repr = (opt->max_iwidth < max_repr) ? opt->max_iwidth : max_repr;
  
  while(start + min_seed_length - 1 < query_len) {
      //fprintf(stderr, "start=%d min_seed_length=%d max_seed_length=%d\n", start, min_seed_length, max_seed_length);
      
      // init iter
      tmap_map4_aux_smem_iter_set_query(iter, max_seed_length, query + start);

      // iterate
      while (0 < tmap_map4_aux_smem_iter_next(iter, bwt)) {
          // go through matches
          for (i = 0; i < iter->matches->n; ++i) {
              tmap_bwt_smem_intv_t *p = &iter->matches->a[i];

              if(p->size <= 0) continue;

              //fprintf(stderr, "EM\t%d\t%d\t%ld\t%llu\t%llu\n", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->size, p->x[0], p->x[1]);
              // too short
              if ((uint32_t)(p->info & 0xFFFF) - (p->info>>32) < min_seed_length) continue;
              // update total
              total++;
              // too many hits?
              if (p->size <= opt->max_iwidth || p->size < max_repr) {
                  // OK
                  tmap_bwt_smem_intv_vec_push(matches, p);
              }
              else if (0 < max_repr) {
                  tmap_bwt_smem_intv_t q;
                  double pr = 0.0;
                  int32_t c;
                  // Keep only representative hits

                  /*
                  p->size = (opt->max_iwidth < max_repr) ? opt->max_iwidth : max_repr;
                  tmap_bwt_smem_intv_vec_push(matches, p);
                  */

                  if (1 == opt->rand_repr) {
                      // choose randomly the representitive hits
                      pr = max_repr / (double)p->size;
                      q = *p;
                      for (k = c = 0; k < q.size; ++k) {
                          if (tmap_rand_get(rand) < pr) {
                              // fake
                              p->x[0] = q.x[0] + k;
                              p->x[1] = q.x[1] + k;
                              p->size = 1;
                              // push
                              tmap_bwt_smem_intv_vec_push(matches, p);
                              // update count
                              c++;
                              if(max_repr <= c) break;
                          }
                      }
                      // reset
                      *p = q;
                  }
                  else {
                      // choose uniformly the representitive hits
                      pr = p->size / (double)max_repr;
                      q = *p;
                      for (k = c = 0; k < q.size; ++k, ++c) {
                          if (pr < c) {
                              // fake
                              p->x[0] = q.x[0] + k;
                              p->x[1] = q.x[1] + k;
                              p->size = 1;
                              // push
                              tmap_bwt_smem_intv_vec_push(matches, p);
                              // update count
                              c = 0;
                          }
                      }
                      // reset
                      *p = q;
                  }
              }
          }
      }
      start += by;
  }
  
  // HERE
  int32_t local_extend = 1;

  // remove seeds if there were too many repetitive hits
  //fprintf(stderr, "matches->n=%d total=%d\n", (int)matches->n, total);
  if (0 < matches->n && !((matches->n / (double)total) < opt->hit_frac)) {
      // realloc
      for (i = n = 0; i < matches->n; ++i) {
          tmap_bwt_smem_intv_t *p = &matches->a[i];
          //fprintf(stderr, "i=%d matches->n=%d p->size=%u\n", i, matches->n, p->size);
          n += p->size;
      }
      tmap_map_sams_realloc(sams, n);
      // go through the matches
      for (i = n = 0; i < matches->n; ++i) {
          tmap_bwt_smem_intv_t *p = &matches->a[i];
          for (k = 0; k < p->size; ++k) {
              tmap_bwt_int_t pacpos;
              tmap_bwt_int_t tmp;
              uint32_t seqid, pos;
              uint8_t strand;
              uint32_t qstart, qend, len;

              qstart = (uint32_t)((p->info >> 32) & 0xFFFF); // should be zero if moved all the way to the start of the query
              qend = (uint32_t)(p->info & 0xFFFF) - 1; // should be qlen if moved all the way to the end of the query
              len = qend - qstart + 1; 
              
              pacpos = p->x[0] + k;
              tmp = pacpos;
              if(bwt->seq_len < pacpos) pacpos = bwt->seq_len;

              // get the packed position
              pacpos = tmap_sa_pac_pos_hash(sa, bwt, pacpos, hash) + 1; // make zero based
              //fprintf(stderr, "X0 p->x[0]=%llu p->x[1]=%llu p->size=%llu k=%llu bwt->seq_len=%llu pacpos=%llu tmp=%llu\n", p->x[0], p->x[1], p->size, k, bwt->seq_len, pacpos, tmp);
              
              // convert to reference co-ordinates
              //if(pacpos == 0 || bwt->seq_len < pacpos) tmap_bug();
              if(0 < tmap_refseq_pac2real(refseq, pacpos, 1, &seqid, &pos, &strand)) {
                  tmap_map_sam_t *s;

                  //fprintf(stderr, "1 seqid:%u pos:%u strand:%d qstart=%u qend=%u len=%u\n", seqid, pos, strand, qstart, qend, len);
                  // adjust the position
                  if(1 == strand) {
                      // NB: could cross a contig boundary 
                      if(pos < len) pos = 0;
                      else pos -= len-1;
                  }
                  //fprintf(stderr, "2 seqid:%u pos:%u strand:%d qstart=%u qend=%u len=%u\n", seqid, pos, strand, qstart, qend, len);

                  // save
                  s = &sams->sams[n];

                  // save the hit
                  s->algo_id = TMAP_MAP_ALGO_MAP4;
                  s->algo_stage = opt->algo_stage;
                  s->strand = strand;
                  s->seqid = seqid;
                  s->pos = pos;
                  s->target_len = query_len;
                  if(refseq->annos[seqid].len < s->target_len) {
                      s->target_len = refseq->annos[seqid].len;
                  }
                  s->score_subo = INT32_MIN;

                  // map4 aux data
                  tmap_map_sam_malloc_aux(s);

                  // TODO forward strand only
                  if(1 == local_extend) {
                      int32_t j, tmp;
                      int32_t score = len * opt->score_match;
                      int32_t score_left, score_right;
                      uint8_t *q_local = NULL, *t_local = NULL;
                      uint32_t q_len_local,t_len_local;
                      uint32_t start_pos, end_pos;
                      int32_t matrix[25];
                      tmap_sw_param_t ap;

                      // init params
                      ap.matrix = matrix;
                      __map_util_gen_ap(ap, opt);
                      
                      score_left = score_right = 0;

                      /*
                      fprintf(stderr, "s->pos=%d score=%d qstart=%d qend=%d len=%d query_len=%d\n", 
                              s->pos, score, qstart, qend, len, query_len);
                              */

                      // extend query left
                      if(0 < qstart) {
                          // query
                          q_local = query;
                          q_len_local = qstart;

                          if(pos <= 1) {
                              // TODO no more bases to the left
                              tmap_bug();
                          }

                          // add in band width
                          t_len_local = q_len_local + opt->bw;
                          if(2 * q_len_local < t_len_local) t_len_local = 2 * q_len_local;

                          // start/end position (one-based)
                          if(0 == strand) {
                              end_pos = pos-1;
                              if(end_pos <= t_len_local) start_pos = 1;
                              else start_pos = end_pos + 1 - t_len_local;
                          }
                          else {
                              start_pos = pos+len;
                              if(refseq->annos[s->seqid].len <= start_pos + t_len_local - 1) end_pos = refseq->annos[s->seqid].len;
                              else end_pos = start_pos + t_len_local - 1;
                          }
                          t_len_local = end_pos - start_pos + 1;

                          // allocate memory
                          t_local = tmap_malloc(sizeof(uint8_t) * t_len_local, "t_local"); // TODO: pre-allocate

                          //fprintf(stderr, "s->seqid=%d start_pos=%d end_pos=%d\n", s->seqid+1, start_pos, end_pos);

                          // NB: IUPAC codes are turned into mismatches
                          if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, t_local, 1, NULL)) {
                              tmap_error("bug encountered", Exit, OutOfRange);
                          }

                          if(1 == strand) { // reverse compliment
                              for(j=0;j<t_len_local>>1;j++) {
                                  tmp = t_local[t_len_local-j-1];
                                  t_local[t_len_local-j-1] = (4 <= t_local[j]) ? t_local[j] : (3-t_local[j]);
                                  t_local[j] = (4 <= tmp) ? tmp : (3-tmp);
                              }
                              if(0 != (t_len_local&1)) t_local[j] = (4 <= t_local[j]) ? t_local[j] : (3-t_local[j]);
                          }

                          score_left = tmap_sw_clipping_core(t_local, t_len_local,
                                                             q_local, q_len_local,
                                                             &ap,
                                                             0, 0, // TODO set these based on strand
                                                             NULL, NULL, // TODO set alignment bounding here
                                                             0); // TODO right justified?
                          
                          //fprintf(stderr, "extend left q_len_local=%d t_len_local=%d score_left=%d\n", q_len_local, t_len_local, score_left);

                          free(t_local);
                          t_local = NULL;
                      }

                      // extend right
                      if(qend < query_len -1) {
                          // query
                          q_local = query + len;
                          q_len_local = query_len - qend - 1;

                          if(refseq->annos[s->seqid].len <= pos) {
                              // TODO no more bases to the right
                              tmap_bug();
                          }

                          // add in band width
                          t_len_local = q_len_local + opt->bw;
                          if(2 * q_len_local < t_len_local) t_len_local = 2 * q_len_local;

                          // start/end position (one-based)
                          if(0 == strand) {
                              start_pos = pos+len;
                              if(refseq->annos[s->seqid].len <= start_pos + t_len_local - 1) end_pos = refseq->annos[s->seqid].len;
                              else end_pos = start_pos + t_len_local - 1;
                          }
                          else {
                              end_pos = pos-1;
                              if(end_pos <= t_len_local) start_pos = 1;
                              else start_pos = end_pos + 1 - t_len_local;
                          }
                          t_len_local = end_pos - start_pos + 1;

                          // allocate memory
                          t_local = tmap_malloc(sizeof(uint8_t) * t_len_local, "t_local"); // TODO: pre-allocate

                          //fprintf(stderr, "s->seqid=%d start_pos=%d end_pos=%d\n", s->seqid+1, start_pos, end_pos);

                          // NB: IUPAC codes are turned into mismatches
                          if(NULL == tmap_refseq_subseq2(refseq, s->seqid+1, start_pos, end_pos, t_local, 1, NULL)) {
                              tmap_error("bug encountered", Exit, OutOfRange);
                          }
                          
                          if(1 == strand) { // reverse compliment
                              for(j=0;j<t_len_local>>1;j++) {
                                  tmp = t_local[t_len_local-j-1];
                                  t_local[t_len_local-j-1] = (4 <= t_local[j]) ? t_local[j] : (3-t_local[j]);
                                  t_local[j] = (4 <= tmp) ? tmp : (3-tmp);
                              }
                              if(0 != (t_len_local&1)) t_local[j] = (4 <= t_local[j]) ? t_local[j] : (3-t_local[j]);
                          }

                          score_right = tmap_sw_clipping_core(t_local, t_len_local,
                                                             q_local, q_len_local,
                                                             &ap,
                                                             0, 0, // TODO set these based on strand
                                                             NULL, NULL, // TODO set alignment bounding here
                                                             0); // TODO right justified?
                          
                          //fprintf(stderr, "extend right q_len_local=%d t_len_local=%d score_right=%d\n", q_len_local, t_len_local, score_right);

                          free(t_local);
                          t_local = NULL;
                      }
                      /*
                      fprintf(stderr, "score=%d score_left=%d score_right=%d\n",
                              score, score_left, score_right);
                              */
                  }

                  n++;
              }
          }
      }
      // realloc
      tmap_map_sams_realloc(sams, n);
  }

  tmap_bwt_smem_intv_vec_destroy(matches);


  return sams;
}
