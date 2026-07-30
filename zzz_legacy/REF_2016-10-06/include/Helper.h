#ifndef HELPER_H
#define HELPER_H

#include <sstream>
#include <string>
#include <vector>
#include <gmpxx.h>

typedef mpq_class Rational;

std::vector<std::string> splitString(const std::string &s, char delim);
size_t getSupportSize(uint64_t support);
size_t getPositionOfLowestSetBit(uint64_t bitsetm);





#endif // HELPER_H
