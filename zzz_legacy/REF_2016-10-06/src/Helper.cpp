#include <Helper.h>

std::vector<std::string> splitString(const std::string &s, char delim) {
    std::stringstream ss;
    std::vector<std::string> elems;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

size_t getSupportSize(uint64_t support)
{
    return __builtin_popcount((unsigned long long)support);
}

size_t getPositionOfLowestSetBit(uint64_t bitsetm) //zero based!!!!
{
    size_t m=0;
    while (true) {
        if ((bitsetm & (1ull << m)) != 0)
            break;
        m++;
    }
    return m;
}








