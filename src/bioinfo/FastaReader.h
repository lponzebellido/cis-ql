#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <string>
#include <vector>

struct FastaRecord {
  std::string header;
  std::string sequenceId;
  std::string sequence;
};

class FastaReader {
public:
  static std::vector<FastaRecord> read(const std::string& filename);
};

#endif
