#ifndef FMAP_REFSEQ_H_
#define FMAP_REFSEQ_H_

#include <stdint.h>
#include "fmap_io.h"

// version id of the indexes
#define FMAP_REFSEQ_VERSION_ID 0
// seed for our random number generator
#define FMAP_REFSEQ_SEED 13
// we use 2 bits per DNA base
#define FMAP_REFSEQ_BYTE_LENGTH(_len) (((_len) + 1) >> 2)
#define FMAP_REFSEQ_BYTE_POSITION(_i) ((3 - ((_i) & 3)) << 1)
// stores the _ith integer base 
#define FMAP_REFSEQ_STORE(_refseq, _i, _c) (_refseq->seq[FMAP_REFSEQ_BYTE_LENGTH(_i)] |= _c << FMAP_REFSEQ_BYTE_POSITION(_i))
// get the _ith integer base
#define FMAP_REFSEQ_BASE(_refseq, _i) ((_refseq->seq[FMAP_REFSEQ_BYTE_LENGTH(_i)] >> FMAP_REFSEQ_BYTE_POSITION(_i)) & 3)
// refseq file extension
#define FMAP_REFSEQ_FILE_EXTENSION ".fmap.refseq"
// compression for the refseq file
#define FMAP_REFSEQ_COMPRESSION FMAP_FILE_NO_COMPRESSION

/*! @typedef
  @abstract  
  @field  name    the name of the contig
  @field  len     the length of the current contig 
  @field  offset  the offset from the start of the reference
  */
typedef struct {
    char *name;
    uint64_t len;
    uint64_t offset;
} fmap_anno_t;

/*! @typedef
  @abstract  
  @field  version_id    the version id of this file
  @field  seed          the random base generator seed
  @field  seq           the packed nucleotide sequence, with contigs concatenated
  @field  annos         the annotations about the contigs
  @field  len           the total length of the reference sequence
  */
typedef struct {
    uint64_t version_id;
    uint32_t seed;
    char *seq;
    fmap_anno_t *annos;
    uint32_t num_annos;
    uint64_t len;
} fmap_refseq_t;

/*! @function
  @abstract
  @field  fn_fasta     the file name of the fasta file
  @field  compression  the type of compression, if any to be used
  @return             a pointer to the initialized memory
 */
fmap_refseq_t *
fmap_refseq_read_fasta(const char *fn_fasta, int32_t compression);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which to store the data 
  @field  prefix   the prefix of the file to be written, usually the fasta file name 
 */
void
fmap_refseq_write(fmap_refseq_t *refseq, const char *prefix);

/*! @function
  @abstract
  @field  prefix   the prefix of the file to be read, usually the fasta file name 
  @return        a pointer to the initialized memory
  */
fmap_refseq_t *
fmap_refseq_read(const char *prefix);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which the data is stored
  */
void
fmap_refseq_destroy(fmap_refseq_t *refseq);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which the data is stored
  @return        the number of bytes used by this data structure
  */
size_t
fmap_refseq_size(fmap_refseq_t *refseq);
#endif
