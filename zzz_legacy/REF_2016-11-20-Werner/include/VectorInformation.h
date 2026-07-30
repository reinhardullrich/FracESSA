#ifndef VECTORINFORMATION_H
#define VECTORINFORMATION_H

#include <vector>
#include <iostream>
#include <string>

#include "Helper.h"

class VectorInformation
{
    public:
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


        std::string getVectorInformation()
        {
            std::string str ="";
            str+=std::to_string(vectorID) + ";";
            for (auto x : resultVector)
                str+= x.template convert_to<std::string>()+",";
            str.pop_back();
            str+= ";" + std::to_string(support) + ";";
            str+= std::to_string(supportSize) + ";";
            str+= std::to_string(extendedSupport) + ";";
            str+= std::to_string(extendedSupportSize) + ";";
            str+= std::to_string(shiftReference) + ";";
            str+= std::to_string(IsEss) + ";";
            str+= std::to_string(ReasonEss) + ";";
            str+= payoff.template convert_to<std::string>() + ";";
            str+= std::to_string(payoffDouble);

            return str;
        }

        static std::string getHeader()
        {
            std::string str ="";
            str+= "VectorID;";
            str+= "Vector;";
            str+= "Support;";
            str+= "SupportSize;";
            str+= "ExtendedSupport;";
            str+= "ExtendedSupportSize;";
            str+= "ShiftReference;";
            str+= "IsEss;";
            str+= "Reason;";
            str+= "Payoff;";
            str+= "PayoffDecimal";

            return str;
        }
};

#endif // VECTORINFORMATION_H
