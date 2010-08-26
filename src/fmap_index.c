#include <stdlib.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"

typedef struct {
    char *fn_fasta;
} fmap_index_opt_t;

void fmap_index_core(fmap_index_opt_t *opt)
{
  fmap_refseq_t *refseq = NULL;

  refseq = fmap_refseq_read_fasta(opt->fn_fasta, FMAP_REFSEQ_COMPRESSION); // read in
  fmap_refseq_write(refseq, opt->fn_fasta); // write out

  fmap_refseq_destroy(refseq); // destroy
}

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s index [optionsn", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "         -f FILE     the FASTA file name\n");
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "\n");
  return 1;
}

int fmap_index(int argc, char *argv[])
{
  int c;
  fmap_index_opt_t opt;

  opt.fn_fasta = NULL;

  while((c = getopt(argc, argv, "f:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'h':
        default:
          return usage();
      }
  }

  if(NULL == opt.fn_fasta) {
      fmap_error("required option -f", Exit, CommandLineArgument);
  }
  if(argc != optind) {
      return usage();
  }

  fmap_index_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
