using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace EFR
{
    class FindEqDouble
    {

        private int _dimension;
        private double[,] _spielMatrixDouble;

        private double[] _outputVector;
        private double[,] _linEqMatrixDouble;
        private double[] _vectorFromLinearEquations;
        private int _IofP;
        private int _supportDimension;
        private int _lastSupportDimension = 1;
        private double _errorBoundSingular = 1e-14; //für singuläre matrizen
        double _errorBoundGauss = 1e-5; //ergebnis des gauss alg.
        double _errorBoundRowSum;

        public FindEqDouble(double[,] spielMatrixDouble)
        {
            _spielMatrixDouble = spielMatrixDouble;
            _dimension = spielMatrixDouble.GetLength(0);
            _outputVector = new double[_dimension];
            _linEqMatrixDouble = new double[_lastSupportDimension + 1, _lastSupportDimension + 2];
            _vectorFromLinearEquations= new double[_lastSupportDimension + 1];
            _errorBoundRowSum = 5 * _dimension * _errorBoundGauss; //laut fehlerfortpflanzung e=e1+e1 sollte es 2n mal fehler von oben sein, 5 als schutz!
        }

        public bool IsFindEq(int IofP, int supportDimension)
        {
            _IofP = IofP;
            _supportDimension = supportDimension;
            if (_supportDimension != _lastSupportDimension)
            {
                _linEqMatrixDouble = new double[_supportDimension + 1, _supportDimension + 2];
                _vectorFromLinearEquations = new double[_supportDimension + 1];
            }
            MakeLinEqMatrix();

            if (!LinEqSolve())
                return false;

            int tracker = 0;
            for (int i = 0; i < _dimension; i++)
            {
                if ((IofP & (1 << i)) != 0)
                {
                    double x = _vectorFromLinearEquations[tracker];
                    if (x > -_errorBoundGauss)
                        _outputVector[i] = x;
                    else
                        return false;
                    tracker += 1;
                }
                else
                    _outputVector[i] = 0D;
            }

            for (int i = 0; i < _dimension; i++)
            {
                if ((IofP & (1 << i)) == 0) //nicht im support iofp - zeilen
                {
                    double rowSum = 0D;
                    for (int j = 0; j < _dimension; j++)
                        if ((IofP & (1 << j)) != 0) // im support von iofp - spalten
                            rowSum += _spielMatrixDouble[i, j] * _outputVector[j];

                    if (!(rowSum <= _vectorFromLinearEquations[supportDimension] + _errorBoundRowSum))
                        return false;
                }

            }
            return true;
        }

        private void MakeLinEqMatrix()
        {

            int row = 0;
            int column = 0;
            for (int i = 0; i < _dimension; i++)
                if ((_IofP & (1 << i)) != 0)
                {
                    column = 0;
                    for (int j = 0; j < _dimension; j++)
                        if ((_IofP & (1 << j)) != 0)
                        {
                            _linEqMatrixDouble[row, column] = _spielMatrixDouble[i, j];
                            column += 1;
                        }
                    _linEqMatrixDouble[row, column] = -1D;
                    _linEqMatrixDouble[row, column + 1] = 0D;
                    row += 1;

                }
            for (int i = 0; i < _supportDimension; i++)
                _linEqMatrixDouble[_supportDimension, i] = 1D;
            _linEqMatrixDouble[_supportDimension, _supportDimension] = 0D;
            _linEqMatrixDouble[_supportDimension, _supportDimension + 1] = 1D;
        }



        private  bool LinEqSolve( )
        {
            int i, j, k, maxrow;
            double tmp, tmp2;
            int n = _supportDimension + 1;

            for (i = 0; i < n; i++)
            {
                /* Find the row with the largest first value */
                maxrow = i;
                for (j = i + 1; j < n; j++)
                    if (Math.Abs(_linEqMatrixDouble[j, i]) > Math.Abs(_linEqMatrixDouble[maxrow, i]))
                        maxrow = j;

                /* Swap the maxrow and ith row */
                for (k = i; k < n + 1; k++)
                {
                    tmp = _linEqMatrixDouble[i, k];
                    _linEqMatrixDouble[i, k] = _linEqMatrixDouble[maxrow, k];
                    _linEqMatrixDouble[maxrow, k] = tmp;
                }

                /* Singular matrix? */
                if (Math.Abs(_linEqMatrixDouble[i, i]) < _errorBoundSingular)
                //if (M[i, i] == 0D)
                {
                    //Debug.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Singuläre Matrix");
                    return false;
                }
                tmp2 = 1D / _linEqMatrixDouble[i, i];
                /* Eliminate the kth element of the jth row */
                for (j = i + 1; j < n; j++)
                {
                    tmp = _linEqMatrixDouble[j, i] * tmp2;
                    for (k = n; k >= i; k--)
                        _linEqMatrixDouble[j, k] -= _linEqMatrixDouble[i, k] * tmp;
                }
            }

            /* Do the back substitution */
            for (j = n - 1; j >= 0; j--)
            {
                tmp = 0;
                for (k = j + 1; k < n; k++)
                    tmp += _linEqMatrixDouble[j, k] * _vectorFromLinearEquations[k];
                _vectorFromLinearEquations[j] = (_linEqMatrixDouble[j, n] - tmp) / _linEqMatrixDouble[j, j];
            }

            return true;
        }

    }
}





//double rowSum = 0D;
//for (int j = 0; j < _spielMatrix.n; j++)
//    rowSum += _spielMatrixDouble[shiftDim, j] * vectorArrayDouble[j];
//if (!(rowSum <= myxDouble[supportDimension] + 0.1D))
//    return false;