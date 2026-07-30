using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
//using System.Threading.Tasks;

namespace EFR
{
    class FindEqBigRational
    {
        public int Support=-1;
        public int SupportSize=-1;
        public int ExtendedSupport=-1;
        public int ExtendedSupportSize = -1;
        public BigRational Payoff;
        public MatrixBigRational ResultVector;

        private MatrixBigRational _spielMatrix;
        private int _dimension;       
        private BigRational[,] _linearEquationMatrix; 
        private BigRational[] _vectorFromLinearEquation;
        private int _lastSupportSize = 1;
        private BigRational _zero;

        public FindEqBigRational(MatrixBigRational spielMatrix)
        {
            _spielMatrix = spielMatrix;
            _dimension = spielMatrix.Rows;
            ResultVector = new MatrixBigRational(_dimension,1);
            _linearEquationMatrix = new BigRational [_lastSupportSize + 1, _lastSupportSize + 2];
            _vectorFromLinearEquation = new BigRational[_lastSupportSize + 1];
            _zero = BigRational.Zero;
            
        }

        public bool FindEquilibrium(int support)
        {
            Support = support;
            SupportSize = BitVector.SupportSize(Support);
            if (SupportSize != _lastSupportSize)
            {
                _linearEquationMatrix = new BigRational[SupportSize + 1, SupportSize + 2];
                _vectorFromLinearEquation = new BigRational[SupportSize + 1];
                _lastSupportSize = SupportSize;
            }
            
            BuildMatrix();

            if (!SolveLinearEquation())
                    return false;

            int tracker = 0;
            for (int i = 0; i < _dimension; i++)
            {
                if ((Support & (1 << i)) != 0)
                {
                    BigRational x = _vectorFromLinearEquation[tracker];
                    if (x.Numerator.Sign == 1)
                        ResultVector[i,0] = x;
                    else
                        return false;
                    tracker += 1;
                }
                else
                    ResultVector[i,0] = _zero;
            }

            ExtendedSupport = Support;
            for (int i = 0; i < _dimension; i++)
            {
                if ((Support & (1 << i)) == 0) //nicht im support iofp - zeilen
                {
                    BigRational rowSum = _zero;
                    for (int j = 0; j < _dimension; j++)
                        if ((Support & (1 << j)) != 0) // im support von iofp - spalten
                            rowSum += _spielMatrix[i, j] * ResultVector[j,0];
                    if (rowSum > _vectorFromLinearEquation[SupportSize])
                        return false;
                    if (rowSum == _vectorFromLinearEquation[SupportSize])
                        ExtendedSupport = (ExtendedSupport | (1 << i));
                }
            }
            ExtendedSupportSize = BitVector.SupportSize(ExtendedSupport);
            Payoff = _vectorFromLinearEquation[SupportSize];
            return true;
        }

        private void BuildMatrix()
        {
            int row = 0;
            int column = 0;
            for (int i = 0; i < _dimension; i++)
                if ((Support & (1 << i)) != 0)
                {
                    column = 0;
                    for (int j = 0; j < _dimension; j++)
                        if ((Support & (1 << j)) != 0)
                        {
                            _linearEquationMatrix[row, column] = _spielMatrix[i, j];
                            column += 1;
                        }
                    _linearEquationMatrix[row, column] = BigRational.MinusOne;
                    _linearEquationMatrix[row, column + 1] = _zero;
                    row += 1;

                }
            for (int i = 0; i <SupportSize; i++)
                _linearEquationMatrix[SupportSize, i] = BigRational.One;
            _linearEquationMatrix[SupportSize, SupportSize] = _zero;
            _linearEquationMatrix[SupportSize, SupportSize + 1] = BigRational.One;
        }



        private bool SolveLinearEquation()
        {
            int i, j, k, maxrow;
            BigRational tmp, tmp2;
            int n = SupportSize + 1;

            for (i = 0; i < n; i++)
            {
                /* Find the row with the largest first value */
                maxrow = i;
                for (j = i + 1; j < n; j++)
                    if (BigRational.Abs(_linearEquationMatrix[j, i]) > BigRational.Abs(_linearEquationMatrix[maxrow, i]))
                        maxrow = j;

                /* Swap the maxrow and ith row */
                for (k = i; k < n + 1; k++)
                {
                    tmp = _linearEquationMatrix[i, k];
                    _linearEquationMatrix[i, k] = _linearEquationMatrix[maxrow, k];
                    _linearEquationMatrix[maxrow, k] = tmp;
                }

                /* Singular matrix? */
                if (_linearEquationMatrix[i, i] == _zero)
                {
                    //Console.WriteLine("!!!!! singular!!!!");
                    return false;
                }

                /* Eliminate the ith element of the jth row */
                tmp2 = BigRational.One / _linearEquationMatrix[i, i];
                for (j = i + 1; j < n; j++)
                {
                    tmp = _linearEquationMatrix[j, i] * tmp2;
                    for (k = n; k >= i; k--)
                        _linearEquationMatrix[j, k] -= _linearEquationMatrix[i, k] * tmp;
                }
            }

            /* Do the back substitution */
            for (j = n - 1; j >= 0; j--)
            {
                tmp = _zero;
                for (k = j + 1; k < n; k++)
                    tmp += _linearEquationMatrix[j, k] * _vectorFromLinearEquation[k];
                _vectorFromLinearEquation[j] = (_linearEquationMatrix[j, n] - tmp) / _linearEquationMatrix[j, j];
            }

            return true;
        }
    }
}


//    Vector Ax = _spielMatrix.MultiplyWithVector(returnVector);
//    for (int shiftDim = 0; shiftDim < _spielMatrix.n; shiftDim++)
//    {
//        if (Ax.GetValue(shiftDim) == v)
//        {
//            extSupport = (extSupport | (1 << shiftDim));
//            extSupportDimension += 1;
//        }
//        if (Ax.GetValue(shiftDim) > v)
//            return null;
//    }
//    return returnVector;
//}