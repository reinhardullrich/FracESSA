#ifndef ESSFINDER_H
#define ESSFINDER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <bitset>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/common_factor.hpp>

#include "MatrixTemplate.h"
#include "VectorInformation.h"
#include "config.h"


class EssFinder
{
    public:
        EssFinder(MatrixTemplate<Rational> matrix);

        void computeEss();
        ~EssFinder();



    private:

        MatrixTemplate<Rational> mGamematrix;
        MatrixTemplate<double> mGamematrixDouble;

        size_t mNumberOfEss = 0;
		size_t mNumberOfCandidates = 0;
		std::ofstream* mLogfile;
		size_t mDimension;
		std::vector< std::vector<uint64_t>> mAllSupports;
		std::vector<bool> mCoprime;
		bool mIsCyclicallySymmetric=false;
		std::vector<VectorInformation> mVectorList;


		size_t searchOneSupportsize(size_t supportSize);
        bool hasCandidateDouble(uint64_t support, size_t supportSize, MatrixTemplate<double> &linearEquationMatrix);
        bool getCandidate(uint64_t support, size_t supportSize, VectorInformation& candidate, MatrixTemplate<Rational> &linearEquationMatrix);
        bool checkStability(uint64_t support, uint64_t extendedSupport, size_t extendedSupportSize, ReasonEss& reasonEss);
};

#endif // ESSFINDER_H


