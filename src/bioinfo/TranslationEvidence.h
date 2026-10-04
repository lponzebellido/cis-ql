#ifndef TRANSLATION_EVIDENCE_H
#define TRANSLATION_EVIDENCE_H

#include <string>

struct TranslationEvidence {
  bool present = false;
  int geneticCode = 0;
  int frame = 0;
  std::string proteinSequence;
};

#endif
