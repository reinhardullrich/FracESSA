#include "VectorInformation.h"

VectorInformation::VectorInformation()
{
    //ctor
}

VectorInformation::~VectorInformation()
{
    //dtor
}

std::string VectorInformation::getVectorInformation()
{
    std::string str ="";
    str+= "'" + std::to_string(vectorID) + "','";
    for (auto x : resultVector)
        str+= x.get_str()+",";
    str.pop_back();
    str+= "','" + std::to_string(support) + "',";
    str+= "'" + std::to_string(supportSize) + "',";
    str+= "'" + std::to_string(extendedSupport) + "',";
    str+= "'" + std::to_string(extendedSupportSize) + "',";
    str+= "'" + std::to_string(shiftReference) + "',";
    str+= "'" + std::to_string(IsEss) + "',";
    str+= "'" + std::to_string(ReasonEss) + "',";
    str+= "'" + payoff.get_str() + "',";
    str+= "'" + std::to_string(payoffDouble) + "',";

    return str;

}


std::string VectorInformation::getHeader()
{
    std::string str ="";
    str+= "'VectorID',";
    str+= "'Vector',";
    str+= "'Support',";
    str+= "'SupportSize',";
    str+= "'ExtendedSupport',";
    str+= "'ExtendedSupportSize',";
    str+= "'ShiftReference',";
    str+= "'IsEss',";
    str+= "'Reason',";
    str+= "'Payoff',";
    str+= "'PayoffDecimal'";

    return str;

}
