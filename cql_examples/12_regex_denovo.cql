
LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS tata_like_anchors;

FIND MOTIF "ATG(?:(?!TAA|TAG|TGA)[ACGT]{3})*(?:TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM tata_like_anchors
  STRAND POSITIVE
  AS in_frame_start_stop_candidates
  WHERE LENGTH >= 15 BP AND LENGTH MOD 3 = 0;

EXPORT in_frame_start_stop_candidates
  TO "in_frame_start_stop_candidates.tsv" FORMAT TSV;
