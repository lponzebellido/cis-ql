#ifndef BED_READER_H
#define BED_READER_H

#include "GenomicRegion.h"
#include <string>
#include <vector>

class BEDReader {
public:
  static std::vector<GenomicRegion> read(const std::string &filename,
                                         const std::string &format,
                                         const std::string &trackAlias,
                                         std::string *error = nullptr);
};

#endif
