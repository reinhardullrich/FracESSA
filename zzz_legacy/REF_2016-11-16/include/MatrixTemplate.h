#ifndef MATRIXTEMPLATE_H
#define MATRIXTEMPLATE_H

#include <iostream>
#include <type_traits>
#include <string>
#include <vector>
#include <cmath>

#include "Helper.h"

template <typename T>
class MatrixTemplate
{

    private:
        std::vector<T> mMatrix;
        size_t mRows;
        size_t mCols;
        size_t mLength;

    public:

        size_t getRows() {return mRows;}
        size_t getCols() {return mCols;}

        void streamMatrix(std::ostream& stream)
        {
            for (size_t i=0;i<mRows;i++) {
                for (size_t j=0;j<mCols;j++)
                    stream<<(*this)(i,j)<<",";
                stream<<std::endl;
            }
        }

        void printMatrix()
        {
            for (size_t i=0;i<mRows;i++) {
                for (size_t j=0;j<mCols;j++)
                    std::cout<<(*this)(i,j)<<",";
                std::cout<<std::endl;
            }
        }

        MatrixTemplate() {}

        MatrixTemplate(size_t rows, int columns)         // Matrix Class constructor
        {
            mRows = rows;
            mCols = columns;
            mLength=rows*columns;
            mMatrix.resize(mLength);
        }

        MatrixTemplate(size_t rows, int columns, T defaultValue)         // Matrix Class constructor
        {
            mRows = rows;
            mCols = columns;
            mLength=rows*columns;
            mMatrix.resize(mLength);
            for (int i =0;i<mLength;i++)
                mMatrix[i]=defaultValue;
        }

        MatrixTemplate(std::vector<T> matrix)
        {
            size_t dimension=std::sqrt(matrix.size());
            mRows = dimension;
            mCols = dimension;
            mLength=dimension*dimension;
            mMatrix.resize(mLength);
            mMatrix=matrix;
        }

        MatrixTemplate(std::vector<T> matrixVector, bool cs)
        {
            size_t dimension=matrixVector.size();
            mRows = dimension;
            mCols = dimension;
            mLength=dimension*dimension;
            mMatrix.resize(mLength);

            std::rotate(matrixVector.begin(),matrixVector.end()-1,matrixVector.end());

            for (size_t i=0; i<mRows;i++) {
                for (size_t j=0;j<mCols;j++) {
                    mMatrix[i*mCols + j]=matrixVector[j];
                }
                std::rotate(matrixVector.begin(),matrixVector.end()-1,matrixVector.end());
            }

        }

        MatrixTemplate(size_t dimension, std::vector<T> matrixHalfVector)
        {
            mRows = dimension;
            mCols = dimension;
            mLength=dimension*dimension;
            mMatrix.resize(mLength);

            std::vector<T> matrixVector;
            matrixVector.reserve(dimension);
            matrixVector.insert(matrixVector.end(), matrixHalfVector.begin(), matrixHalfVector.end());

            if (dimension %2 ==0)  //dimension is even
                matrixVector.insert(matrixVector.end(), matrixHalfVector.rbegin()+1, matrixHalfVector.rend());
            else
                matrixVector.insert(matrixVector.end(), matrixHalfVector.rbegin(), matrixHalfVector.rend());

            matrixVector.push_back(0);

            std::rotate(matrixVector.begin(),matrixVector.end()-1,matrixVector.end());

            for (size_t i=0; i<mRows;i++) {
                for (size_t j=0;j<mCols;j++) {
                    mMatrix[i*mCols + j]=matrixVector[j];
                }
                std::rotate(matrixVector.begin(),matrixVector.end()-1,matrixVector.end());
            }

        }

        ~MatrixTemplate()
        {

        }

        size_t rows()
        {
            return mRows;
        }

        size_t cols()
        {
            return mCols;
        }

        // Access the individual elements
        inline T& operator()(size_t row, size_t col)
        {
          return mMatrix[row*mCols+col];
        }

        // Access the individual elements (const)
        inline const T& operator()(size_t row, size_t col) const
        {
          return mMatrix[row*mCols+col];
        }

        MatrixTemplate<double> toDouble()
        {
            MatrixTemplate<double> A = MatrixTemplate<double> (mRows, mCols);
            for (size_t i=0; i<mRows; i++)
                for (size_t j=0; j<mCols;j++) {
                    A(i,j)=static_cast<double>(mMatrix[i*mCols+j]);
                }
            return A;
        }

        inline static MatrixTemplate<T> initializeLinearEquationMatrix(size_t supportSize)
        {
            MatrixTemplate<T> linearEquationMatrix= MatrixTemplate<T>(supportSize+1, supportSize+2);
            return linearEquationMatrix;
        }

        void getLinearEquationMatrix(uint64_t support, size_t supportSize, MatrixTemplate<T> &linearEquationMatrix)
        {
            size_t row = 0;
            size_t column = 0;

            for (size_t i = 0; i < mRows; i++)
                if ((support & (1ull << i)) != 0)
                {
                    column = 0;
                    for (size_t j = 0; j < mCols; j++)
                        if ((support & (1ull << j)) != 0)
                        {
                            linearEquationMatrix(row, column)= mMatrix[i*mCols+j];
                            column += 1;
                        }
                    linearEquationMatrix(row, column) = -1;
                    linearEquationMatrix(row, column+1) = 0;
                    row += 1;

                }
            for (size_t i = 0; i < supportSize; i++)
                linearEquationMatrix(supportSize, i) = 1;

            linearEquationMatrix(supportSize, supportSize) = 0;
            linearEquationMatrix(supportSize, supportSize+1) = 1;
        }

        static MatrixTemplate<T> ZeroMatrix(int rows, int cols)       // Function generates the zero matrix
        {
            MatrixTemplate<T> matrix = MatrixTemplate<T>(rows, cols);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    matrix(i,j) = 0;
            return matrix;
        }

        static MatrixTemplate<T> IdentityMatrix(size_t dimension)   // Function generates the identity matrix
        {
            MatrixTemplate<T> matrix = ZeroMatrix(dimension, dimension);
            for (size_t i = 0; i < dimension; i++)
                matrix(i,i) = 1;
            return matrix;
        }

        bool IsPosDefDouble() //cholesky-decomposition for double!
        {
            int n = mRows;
            MatrixTemplate<double> L = ZeroMatrix(n, n);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    double aSum = 0.;
                    for (int k = 0; k < j; k++)
                        aSum += L(i,k) * L(j,k);
                    L(i,j) = (1 / L(j,j)) * ((*this)(i,j) - aSum);
                }
                double bSum = 0.;
                for (int k = 0; k < i; k++)
                    bSum += L(i,k) * L(i,k);
                double x = (*this)(i,i) - bSum;
                if (x< 0.01)
                    return false;
                L(i,i) = sqrt(x);
            }
            return true;
        }

//		//from sympy/matrices/dense.py
//		L = zeros(self.rows, self.rows)
//			for i in range(self.rows):
//				for j in range(i):
//					L[i, j] = (1 / L[j, j])*(self[i, j] -sum(L[i, k]*L[j, k] for k in range(j)))
//				L[i, i] = sqrt(self[i, i] -sum(L[i, k]**2 for k in range(i)))
//			return self._new(L)


        bool IsPosDef() //ldlt-decomposition for rational numbers
        {
            int n = mRows;
            MatrixTemplate<Rational> D = ZeroMatrix(n, n);
            MatrixTemplate<Rational> L = IdentityMatrix (n);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    Rational aSum = Rational(0);
                    for (int k = 0; k < j; k++)
                        aSum += L(i,k) * L(j,k) * D(k,k);
                    L(i,j) = (1 / D(j,j)) * ((*this)(i,j) - aSum);
                }
                Rational bSum = Rational(0);
                for (int k = 0; k < i; k++)
                    bSum += L(i,k) * L(i,k) * D(k,k);
                D(i,i) = (*this)(i,i) - bSum;
                if (D(i,i)<=Rational(0))
                    return false;
            }
            return true;
        }

//		//from sympy/matrices/dense.py
//		D = zeros(self.rows, self.rows)
//			L = eye(self.rows)
//			for i in range(self.rows):
//				for j in range(i):
//					L[i, j] = (1 / D[j, j])*(self[i, j] - sum(L[i, k]*L[j, k]*D[k, k] for k in range(j)))
//				D[i, i] = self[i, i] - sum(L[i, k]**2*D[k, k] for k in range(i))
//				return self._new(L), self._new(D)



//        /* Recursive function for finding determinant of matrix.
//           n is current dimension of mat[][]. */
//        T determinantLaplace()
//        {
//            T D = 0;
//
//            if (this->mRows==1)
//                return this->mMatrix[0];
//
//            int sign = 1;  // To store sign multiplier
//
//             // Iterate for each element of first row
//            for (size_t f = 0; f < this->mRows; f++)
//            {
//                // Getting Cofactor of mat[0][f]
//                MatrixTemplate<T> temp = this->developBy(0,f);
//                D += sign * (*this)(0,f) * temp.determinantLaplace();
//                // terms are to be added with alternate sign
//                sign = -sign;
//            }
//
//            return D;
//        }


        MatrixTemplate<T> clone()
        {
            MatrixTemplate<T> A = MatrixTemplate<T> (mRows, mCols);
            for (size_t i=0; i<mRows; i++)
                for (size_t j=0; j<mCols;j++) {
                    A(i,j)=mMatrix[i*mCols+j];
                }
            return A;
        }

        T determinant() //lu-decomposition, crout-algorithm from http://algorithm.narod.ru/ln/crout.c
        {
            MatrixTemplate<T> a= this->clone();
            size_t n = mRows;
            std::vector<size_t> indx=std::vector<size_t>(n);
            std::vector<T> vv=std::vector<T>(n);
            int d =1;

            size_t i,imax,j,k;
            T big,sum,temp;

            /* search for the largest element in each row; save the scaling in the
            temporary array vv and return zero if the matrix is singular */
            for(i=0;i<n;i++) {
                big=0.;
                for(j=0;j<n;j++) if((temp=abs(a(i,j)))>big) big=temp;
                if(big==0) return(0);
                vv[i]=big;
             }
             /* the main loop for the Crout's algorithm */
             for(j=0;j<n;j++) {
                /* this is the part a) of the algorithm except for i==j */
                for(i=0;i<j;i++) {
                    sum=a(i,j);
                    for(k=0;k<i;k++) sum-=a(i,k)*a(k,j);
                    a(i,j)=sum;
                }
                /* initialize for the search for the largest pivot element */
                big=0;imax=j;
                /* this is the part a) for i==j and part b) for i>j + pivot search */
                for(i=j;i<n;i++) {
                    sum=a(i,j);
                    for(k=0;k<j;k++) sum-=a(i,k)*a(k,j);
                    a(i,j)=sum;
                    /* is the figure of merit for the pivot better than the best so far? */
                    if((temp=vv[i]*abs(sum))>=big) {big=temp;imax=i;}
                }
                /* interchange rows, if needed, change parity and the scale factor */
                if(imax!=j) {
                    for(k=0;k<n;k++) {temp=a(imax,k);a(imax,k)=a(j,k);a(j,k)=temp;}
                    d=-d;vv[imax]=vv[j];
                }
                /* store the index */
                indx[j]=imax;
                if(a(j,j)==0) return(0);
                /* finally, divide by the pivot element */
                if(j<n-1) {
                    temp=1/a(j,j);
                    for(i=j+1;i<n;i++) a(i,j)*=temp;
                }
            }
            T res=d;
            for(j=0;j<n;j++) res*=a(j,j);
            return(res);
        }

        MatrixTemplate<T> transpose()
        {
            MatrixTemplate<T> A=MatrixTemplate<T>(mRows,mCols);
            for (size_t i=0; i<mRows; i++)
                for (size_t j=0; j<mCols; j++)
                    A(i,j)=(*this)(j,i);
            return A;
        }


        MatrixTemplate<T> developBy(size_t row, size_t col)
        {
            MatrixTemplate<T> A=MatrixTemplate<T>(mRows-1,mCols-1);
            for (size_t i=0; i<mRows-1; i++)
                for (size_t j=0; j<mCols-1; j++) {
                    size_t rowIndexReal = (i>=row) ? i+1 : i;
                    size_t colIndexReal = (j>=col) ? j+1 : j;
                    A(i,j)=(*this)(rowIndexReal,colIndexReal);
                }
            return A;
        }

        MatrixTemplate<T> adjugate()
        {
            MatrixTemplate<T> A=MatrixTemplate<T>(mRows,mCols);

            if (mRows==1) {
                A(0,0)=1;
                return A;
            }

            for (size_t i=0; i<mRows; i++) {
                for (size_t j=0; j<mCols; j++) {
                    A(i,j)=std::pow(-1,i+j) * (this->developBy(i,j)).determinant();
                }

            }
            return A.transpose();
        }

//        MatrixTemplate<T> adjugateLaplace()
//        {
//            MatrixTemplate<T> A=MatrixTemplate<T>(mRows,mCols);
//
//            if (mRows==1) {
//                A(0,0)=1;
//                return A;
//            }
//
//            for (size_t i=0; i<mRows; i++) {
//                for (size_t j=0; j<mCols; j++) {
//                    A(i,j)=std::pow(-1,i+j) * (this->developBy(i,j)).determinantLaplace(); //determinantLaplace()
//                }
//
//            }
//            return A.transpose();
//        }

        MatrixTemplate<T> multiplyWith(MatrixTemplate<T> B)
        {
            MatrixTemplate<T> C= MatrixTemplate<T> (mRows,B.getCols());
            for (size_t i=0; i<mRows; i++) {
                for (size_t j=0; j<B.getCols(); j++) {
                    T sum=0;
                    for (size_t k=0;k<mCols;k++)
                        sum +=(*this)(i,k)*B(k,j);
                    C(i,j)=sum;
                }
            }
            return C;
        }

        bool greaterZero()
        {
            for (size_t i=0; i<mRows; i++)
                for (size_t j=0; j<mCols; j++)
                    if ((*this)(i,j)<=0)
                        return 0;
            return 1;
        }


        bool isCopositive()
        {
//            std::ofstream *mLogfile= new std::ofstream("cp.txt");
            for (uint64_t support = 1ull; support < (1ull << mRows); support++) { //iterate all subsets without the empty set
                size_t n = getSupportSize(support);
                MatrixTemplate<T> A = MatrixTemplate<T>(n,n);

                size_t row = 0;
                size_t column = 0;
                for (size_t i = 0;i<mRows;i++) {
                    if ((support & (1ull << i)) != 0) {
                        column = 0;
                        for (size_t j = 0; j < mCols; j++)
                            if ((support & (1ull << j)) != 0) {
                                A(row,column) = (*this)(i,j);
                                column += 1;
                            }
                        row += 1;
                    }
                }
//                *mLogfile<<"-------------------------------------------------------------------"<<std::endl;
//                A.streamMatrix(*mLogfile);
//                *mLogfile<<A.determinant()<<std::endl;
//                *mLogfile<<A.adjugate().greaterZero()<<std::endl;
                if (A.determinant()<=0 && A.adjugate().greaterZero()) {
                    return 0;
                }

            }
            return 1;
        }

//        bool isCopositiveLaplace()
//        {
//            std::ofstream *mLogfile= new std::ofstream("cpLaplace.txt");
//            for (uint64_t support = 1ull; support < (1ull << mRows); support++) { //iterate all subsets without the empty set
//                size_t n = getSupportSize(support);
//                MatrixTemplate<T> A = MatrixTemplate<T>(n,n);
//
//                size_t row = 0;
//                size_t column = 0;
//                for (size_t i = 0;i<mRows;i++) {
//                    if ((support & (1ull << i)) != 0) {
//                        column = 0;
//                        for (size_t j = 0; j < mCols; j++)
//                            if ((support & (1ull << j)) != 0) {
//                                A(row,column) = (*this)(i,j);
//                                column += 1;
//                            }
//                        row += 1;
//                    }
//                }
//                *mLogfile<<"-------------------------------------------------------------------"<<std::endl;
//                A.streamMatrix(*mLogfile);
//                *mLogfile<<A.determinantLaplace()<<std::endl;
//                *mLogfile<<A.adjugateLaplace().greaterZero()<<std::endl;
//                if (A.determinantLaplace()<=0 && A.adjugateLaplace().greaterZero()) { //determinantLaplace
//                    return 0;
//                }
//
//            }
//            return 1;
//        }

};

#endif // MATRIXTEMPLATE_H

//template<typename T>
//template<class Q = T>
//typename std::enable_if<std::is_same<Q, double>::value, void>::type MatrixTemplate::initNumbers()
//{
//    mZero=0.;
//    mOne=1.;
//}
//
//template<typename T>
//template<class Q = T>
//typename std::enable_if<std::is_same<Q, Rational>::value, void>::type MatrixTemplate::initNumbers()
//{
//    mZero=Rational(0);
//    mOne=Rational(1);
//}


//    public MatrixBigRational(string[] matrixString)
//    {
//        int nrow=(int)Math.Sqrt(matrixString.Length);
//        if (nrow * nrow != matrixString.Length) {
//            Console.WriteLine ("Error: The given gamematrix is not square!!!");
//            return;
//        }
//
//        Rows = nrow;
//        Columns=matrixString.Length/nrow;
//        _matrix = new BigRational[Rows, Columns];
//
//        for (int i = 0; i < Rows; i++)
//            for (int j = 0; j < Columns; j++)
//                this[i, j] = new BigRational(matrixString[i * nrow + j]);
//    }



//
//
//    public override string ToString()
//    {
//        StringBuilder output = new StringBuilder("");
//        for (int i = 0; i < Rows; i++)
//            for (int j = 0; j < Columns; j++)
//            {
//                if (_matrix[i,j].Denominator==1)
//                    output.Append(_matrix[i, j].Numerator).Append(",");
//                else
//                    output.Append(_matrix[i, j]).Append(",");
//            }
//        output.Length--;
//        return output.ToString();
//    }
//
//    public string ToMatrixRepresenation()
//    {
//        StringBuilder output = new StringBuilder("");
//        for (int i = 0; i < Rows; i++) {
//            for (int j = 0; j < Columns; j++) {
//                if (_matrix [i, j].Denominator == 1)
//                    output.Append (_matrix [i, j].Numerator).Append ("\t");
//                else
//                    output.Append (_matrix [i, j]).Append ("\t");
//            }
//            output.Append ("\r\n");
//        }
//        output.Length--;
//        output.Length--;
//        return output.ToString();
//    }
//
//    public double[,] ToDoubleArray()
//    {
//        double[,] output = new double[Rows, Columns];
//        for (int i = 0; i < Rows; i++)
//            for (int j = 0; j < Columns; j++)
//                output[i, j] = ((double)this[i, j]);
//        return output;
//    }
//
//
//

//
//    #endregion public methods
//
//    #region private static methods
//
//    private static MatrixBigRational Add(MatrixBigRational m1, MatrixBigRational m2)
//    {
//        MatrixBigRational r = new MatrixBigRational(m1.Rows, m1.Columns);
//        for (int i = 0; i < r.Rows; i++)
//            for (int j = 0; j < r.Columns; j++)
//                r[i, j] = m1[i, j] + m2[i, j];
//        return r;
//    }
//
//    public static MatrixBigRational Multiply(MatrixBigRational m1, MatrixBigRational m2)                  //  matrix multiplication
//    {
//        MatrixBigRational result = ZeroMatrix(m1.Rows, m2.Columns);
//        for (int i = 0; i < result.Rows; i++)
//            for (int j = 0; j < result.Columns; j++)
//                for (int k = 0; k < m1.Columns; k++)
//                    result[i, j] += m1[i, k] * m2[k, j];
//        return result;
//    }
//    public static MatrixBigRational Multiply(BigRational alpha, MatrixBigRational m)                          // Multiplication by constant alpha
//    {
//        MatrixBigRational r = new MatrixBigRational(m.Rows, m.Columns);
//        for (int i = 0; i < m.Rows; i++)
//            for (int j = 0; j < m.Columns; j++)
//                r[i, j] = m[i, j] * alpha;
//        return r;
//    }
//
//    #endregion private static methods
//
//    #region operators
//
//    public static MatrixBigRational operator -(MatrixBigRational m)
//    { return MatrixBigRational.Multiply(-1, m); }
//
//    public static MatrixBigRational operator +(MatrixBigRational m1, MatrixBigRational m2)
//    { return MatrixBigRational.Add(m1, m2); }
//
//    public static MatrixBigRational operator -(MatrixBigRational m1, MatrixBigRational m2)
//    { return MatrixBigRational.Add(m1, -m2); }
//
//    public static MatrixBigRational operator *(MatrixBigRational m1, MatrixBigRational m2)
//    { return MatrixBigRational.Multiply(m1, m2); }
//
//    public static MatrixBigRational operator *(BigRational alpha, MatrixBigRational m)
//    { return MatrixBigRational.Multiply(alpha, m); }
//


//
//        template<class Q = T>
//        typename std::enable_if<std::is_same<Q, double>::value, std::string>::type  toString(Rational number)
//        {
//            return std::to_string(number);
//        }
//
//        template<class Q = T>
//        typename std::enable_if<std::is_same<Q, Rational>::value, std::string>::type toString(Rational number)
//        {
//            return number.get_str();
//        }


