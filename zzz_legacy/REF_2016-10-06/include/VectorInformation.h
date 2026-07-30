#ifndef VECTORINFORMATION_H
#define VECTORINFORMATION_H

#include <gmpxx.h>

#include <vector>
#include <iostream>
#include <string>

#include <Helper.h>

enum class ReasonEss : int
{

    indeterminate=-1,
    true_pure_ess=1,
    true_posdef_double = 2,
    true_posdef_rational=3,
    true_copositive=4,
    false_not_posdef_and_K_0_1 = 5,
    false_not_partial_copositive = 6,
    false_not_copositive=7
};

class VectorInformation
{
    public:
        VectorInformation();
        virtual ~VectorInformation();
        std::string getVectorInformation();
        static std::string getHeader();

        size_t vectorID;
        std::vector<Rational> resultVector;
        uint64_t support;
        size_t supportSize;
        uint64_t extendedSupport;
        size_t extendedSupportSize;
        size_t shiftReference;
        bool IsEss;
        int ReasonEss;
        Rational payoff;
        double payoffDouble;



};

#endif // VECTORINFORMATION_H
