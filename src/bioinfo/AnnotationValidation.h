#ifndef ANNOTATION_VALIDATION_H
#define ANNOTATION_VALIDATION_H

#include <cstddef>
#include <string>
#include <vector>

struct AnnotationMultiRecordIdentity {
  std::string id;
  std::vector<size_t> recordNumbers;
};

struct AnnotationUnresolvedParent {
  size_t recordNumber = 0;
  std::string chr;
  std::string type;
  std::string childId;
  std::string childName;
  std::string parentId;
};

struct AnnotationIdentityConflict {
  std::string id;
  std::vector<size_t> recordNumbers;
  std::vector<std::string> chromosomes;
  std::vector<std::string> types;
  std::vector<std::string> strands;
};

struct AnnotationValidationReport {
  std::string annotationAlias;
  size_t totalRecords = 0;
  size_t recordsWithId = 0;
  size_t uniqueIds = 0;
  size_t parentReferences = 0;
  std::vector<AnnotationMultiRecordIdentity> multiRecordIds;
  std::vector<AnnotationUnresolvedParent> unresolvedParents;
  std::vector<AnnotationIdentityConflict> identityConflicts;
  std::vector<std::vector<std::string>> cycles;

  bool valid() const {
    return unresolvedParents.empty() && identityConflicts.empty() &&
           cycles.empty();
  }
};

#endif
