#ifndef ANNOTATION_EVIDENCE_H
#define ANNOTATION_EVIDENCE_H

#include <map>
#include <string>
#include <vector>

struct AnnotationEvidence {
  bool present = false;
  std::string source;
  std::string score;
  std::string phase;
  std::string id;
  std::string name;
  std::vector<std::string> parents;
  std::map<std::string, std::string> attributes;
  std::string rawAttributes;
};

#endif
