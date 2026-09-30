#ifndef ANNOTATION_FEATURE_GROUP_H
#define ANNOTATION_FEATURE_GROUP_H

#include "GenomicRegion.h"
#include <cstddef>
#include <string>
#include <vector>

struct AnnotationFeatureGroup {
  std::string sourceAlias;
  std::string id;
  std::string chr;
  std::string type;
  std::string strand;
  size_t spanStart = 0;
  size_t spanEnd = 0;
  size_t totalLength = 0;
  std::vector<GenomicRegion> members;
};

#endif
