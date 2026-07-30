using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Diagnostics;


namespace EFR
{
    class MatrixBigRational
    {

        public int Rows { get; set; }
        public int Columns { get; set; }

        private BigRational[,] _matrix;

        #region constructors

        public MatrixBigRational(int rows, int columns)         // Matrix Class constructor
        {
            Rows = rows;
            Columns = columns;
            _matrix = new BigRational[Rows, Columns];
        }

        public MatrixBigRational(int dimension)         // Matrix Class constructor
        {
            Rows = dimension;
            Columns = dimension;
            _matrix = new BigRational[Rows, Columns];
        }

        public MatrixBigRational(int dimension, string[] circularLineHalf)
        //circularLineHalf for the short form n/2 or n/2 + 1 !!!
        {
            BigRational[] matrixArray = new BigRational[dimension];
            matrixArray[0] = BigRational.Zero;
            BigRational tempBR;

            if (dimension % 2 != 0) //is odd
            {
                for (int i = 0; i < circularLineHalf.Length; i++)
                {
                    tempBR = new BigRational(circularLineHalf[i]);
                    matrixArray[i + 1] = tempBR;
                    matrixArray[dimension - 1 - i] = tempBR;
                }
            }
            else
            {
                for (int i = 0; i < circularLineHalf.Length - 1; i++)
                {
                    tempBR = new BigRational(circularLineHalf[i]);
                    matrixArray[i + 1] = tempBR;
                    matrixArray[dimension - 1 - i] = tempBR;
                }
                matrixArray[dimension / 2] = new BigRational(circularLineHalf[circularLineHalf.Length - 1]);
            }

            Rows = dimension;
            Columns = dimension;
            _matrix = new BigRational[Rows, Columns];
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                    _matrix[i, j] = matrixArray[j];
                Helper.ShiftRight(matrixArray);
            }

        }

        public MatrixBigRational(string[] matrixString)
        {
            int nrow = (int)Math.Sqrt(matrixString.Length);
            if (nrow * nrow != matrixString.Length)
            {
                Console.WriteLine("Error: The given gamematrix is not square!!!");
                return;
            }

            Rows = nrow;
            Columns = matrixString.Length / nrow;
            _matrix = new BigRational[Rows, Columns];

            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    this[i, j] = new BigRational(matrixString[i * nrow + j]);
        }

        #endregion constructors

        #region static methods

        public static MatrixBigRational ZeroMatrix(int iRows, int iCols)       // Function generates the zero matrix
        {
            MatrixBigRational matrix = new MatrixBigRational(iRows, iCols);
            BigRational brZero = BigRational.Zero;
            for (int i = 0; i < iRows; i++)
                for (int j = 0; j < iCols; j++)
                    matrix[i, j] = brZero;
            return matrix;
        }

        public static MatrixBigRational IdentityMatrix(int n)   // Function generates the identity matrix
        {
            MatrixBigRational matrix = ZeroMatrix(n, n);
            for (int i = 0; i < Math.Min(n, n); i++)
                matrix[i, i] = BigRational.One;
            return matrix;
        }

        #endregion static methods

        #region public methods



        public override string ToString()
        {
            StringBuilder output = new StringBuilder("");
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                {
                    if (_matrix[i, j].Denominator == 1)
                        output.Append(_matrix[i, j].Numerator).Append(",");
                    else
                        output.Append(_matrix[i, j]).Append(",");
                }
            output.Length--;
            return output.ToString();
        }

        public string ToMatrixRepresenation()
        {
            StringBuilder output = new StringBuilder("");
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    if (_matrix[i, j].Denominator == 1)
                        output.Append(_matrix[i, j].Numerator).Append("\t");
                    else
                        output.Append(_matrix[i, j]).Append("\t");
                }
                output.Append("\r\n");
            }
            output.Length--;
            output.Length--;
            return output.ToString();
        }

        public double[,] ToDoubleArray()
        {
            double[,] output = new double[Rows, Columns];
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    output[i, j] = ((double)this[i, j]);
            return output;
        }

        //        public bool IsSquare()
        //        {
        //            return (Rows == Columns);
        //        }

        public BigRational this[int iRow, int iCol]      // Access this matrix as a 2D array
        {
            get { return _matrix[iRow, iCol]; }
            set { _matrix[iRow, iCol] = value; }
        }

        public bool IsPosDefDouble()
        {
            int n = this.Rows;
            double[][] A = new double[n][];

            for (int i = 0; i < n; i++)
            {
                double[] temp = new double[n];
                for (int j = 0; j < n; j++)
                    temp[j] = (double)this[i, j];
                A[i] = temp;
            }

            double[][] L = new double[n][];
            for (int i = 0; i < n; i++)
                L[i] = new double[n];

            // Main loop.
            for (int j = 0; j < n; j++)
            {
                double[] Lrowj = L[j];
                double d = 0.0;
                for (int k = 0; k < j; k++)
                {
                    double[] Lrowk = L[k];
                    double s = 0.0;
                    for (int i = 0; i < k; i++)
                    {
                        s += Lrowk[i] * Lrowj[i];
                    }
                    Lrowj[k] = s = (A[j][k] - s) / L[k][k];
                    d = d + s * s;
                    if (A[k][j] != A[j][k])
                        Console.WriteLine("Matrix ist nicht symmetrisch in Cholesky-Decomp.!!!");
                }
                d = A[j][j] - d;
                if (d <= 0.01D)
                    return false;
                L[j][j] = System.Math.Sqrt(d);
                for (int k = j + 1; k < n; k++)
                    L[j][k] = 0.0;
            }
            return true;
        }

        public bool IsPosDefCholesky_TEST()
        {
            int n = this.Rows;
            MatrixBigRational L = MatrixBigRational.ZeroMatrix(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    BigRational aSum = BigRational.Zero;
                    for (int k = 0; k < j; k++)
                        aSum += L[i, k] * L[j, k];
                    L[i, j] = (1 / L[j, j]) * (this[i, j] - aSum);
                }
                BigRational bSum = BigRational.Zero;
                for (int k = 0; k < i; k++)
                    bSum += L[i, k] * L[i, k];
                BigRational x = this[i, i] - bSum;
                if (x.Numerator.Sign != 1)
                    return false;
                L[i, i] = new BigRational(Math.Sqrt((double)x));
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

        public bool IsPosDef()
        {
            int n = this.Rows;
            MatrixBigRational D = MatrixBigRational.ZeroMatrix(n, n);
            MatrixBigRational L = MatrixBigRational.IdentityMatrix(n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    BigRational aSum = BigRational.Zero;
                    for (int k = 0; k < j; k++)
                        aSum += L[i, k] * L[j, k] * D[k, k];
                    //Console.WriteLine (D[j, j]);
                    L[i, j] = (1 / D[j, j]) * (this[i, j] - aSum);
                }
                BigRational bSum = BigRational.Zero;
                for (int k = 0; k < i; k++)
                    bSum += L[i, k] * L[i, k] * D[k, k];
                D[i, i] = this[i, i] - bSum;
                if (D[i, i].Numerator.Sign != 1)
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


        #endregion public methods

        #region private static methods

        private static MatrixBigRational Add(MatrixBigRational m1, MatrixBigRational m2)
        {
            MatrixBigRational r = new MatrixBigRational(m1.Rows, m1.Columns);
            for (int i = 0; i < r.Rows; i++)
                for (int j = 0; j < r.Columns; j++)
                    r[i, j] = m1[i, j] + m2[i, j];
            return r;
        }

        public static MatrixBigRational Multiply(MatrixBigRational m1, MatrixBigRational m2)                  //  matrix multiplication
        {
            MatrixBigRational result = ZeroMatrix(m1.Rows, m2.Columns);
            for (int i = 0; i < result.Rows; i++)
                for (int j = 0; j < result.Columns; j++)
                    for (int k = 0; k < m1.Columns; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }
        public static MatrixBigRational Multiply(BigRational alpha, MatrixBigRational m)                          // Multiplication by constant alpha
        {
            MatrixBigRational r = new MatrixBigRational(m.Rows, m.Columns);
            for (int i = 0; i < m.Rows; i++)
                for (int j = 0; j < m.Columns; j++)
                    r[i, j] = m[i, j] * alpha;
            return r;
        }

        #endregion private static methods

        #region operators

        public static MatrixBigRational operator -(MatrixBigRational m)
        { return MatrixBigRational.Multiply(-1, m); }

        public static MatrixBigRational operator +(MatrixBigRational m1, MatrixBigRational m2)
        { return MatrixBigRational.Add(m1, m2); }

        public static MatrixBigRational operator -(MatrixBigRational m1, MatrixBigRational m2)
        { return MatrixBigRational.Add(m1, -m2); }

        public static MatrixBigRational operator *(MatrixBigRational m1, MatrixBigRational m2)
        { return MatrixBigRational.Multiply(m1, m2); }

        public static MatrixBigRational operator *(BigRational alpha, MatrixBigRational m)
        { return MatrixBigRational.Multiply(alpha, m); }

        #endregion operators


        public bool IsCopositive()
        {
            LinkedList<int> subsetList = BitVector.GetAllSubsetsWithMoreThanZeroElements(this.Rows);
            foreach (int subset in subsetList)
            {
                int dimension = BitVector.SupportSize(subset);
                MatrixBigRational subsetMatrix = this.GetPrincipalSubmatrixByBitmask(subset, dimension);

                if (dimension == 2)
                {
                    if (!(subsetMatrix[0, 0] > BigRational.Zero && subsetMatrix[1, 1] > BigRational.Zero
                        && ((subsetMatrix[0, 0] * subsetMatrix[1, 1] - subsetMatrix[0, 1] * subsetMatrix[1, 0] > BigRational.Zero) || (subsetMatrix[0, 1] >= BigRational.Zero))))
                        return false;
                }
                else
                {
                    //if (subsetMatrix.DeterminantIsSmallerEqualZero() && subsetMatrix.AdjugateMatrixIsGreaterZero())
                        //return false;

                    if (subsetMatrix.Determinant() <= new BigInteger(0) && subsetMatrix.Adjugate().GreaterZero())
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public MatrixBigRational GetPrincipalSubmatrixByBitmask(int bitmask, int length)
        {
            MatrixBigRational mbr = new MatrixBigRational(length);
            int row = 0;
            int column = 0;
            for (int i = 0; i < Rows; i++)
                if ((bitmask & (1 << i)) != 0)
                {
                    column = 0;
                    for (int j = 0; j < Columns; j++)
                        if ((bitmask & (1 << j)) != 0)
                        {
                            mbr[row, column] = _matrix[i, j];
                            column += 1;
                        }
                    row += 1;
                }
            return mbr;
        }

        private MatrixBigRational DevelopBy(int row, int col)
        {
            MatrixBigRational A = new MatrixBigRational(this.Rows - 1, this.Columns - 1);
            for (int i = 0; i < this.Rows - 1; i++)
                for (int j = 0; j < this.Columns - 1; j++)
                {
                    int rowIndexReal = (i >= row) ? i + 1 : i;
                    int colIndexReal = (j >= col) ? j + 1 : j;
                    A[i, j] = this._matrix[rowIndexReal, colIndexReal];
                }
            return A;
        }

        public BigRational Determinant()
        {
            BigRational D = new BigRational(new BigInteger(0));

            if (this.Rows==1)
                return this._matrix[0,0];

            int sign = 1;  // To store sign multiplier

             // Iterate for each element of first row
            for (int f = 0; f < this.Rows; f++)
            {
                // Getting Cofactor of mat[0][f]
                MatrixBigRational temp = this.DevelopBy(0,f);
                D += sign * this._matrix[0,f] * temp.Determinant();
                // terms are to be added with alternate sign
                sign = -sign;
            }

            return D;
        }

        public MatrixBigRational Adjugate()
        {
            MatrixBigRational A = new MatrixBigRational(this.Rows, this.Columns);

            if (this.Rows == 1)
            {
                A[0, 0] = 1;
                return A;
            }

            for (int i = 0; i < this.Rows; i++)
            {
                for (int j = 0; j < this.Columns; j++)
                {
                    A[i, j] = Math.Pow(-1, i + j) * (this.DevelopBy(i, j)).Determinant();
                }

            }
            return A.Transpose();
        }

        public MatrixBigRational Transpose()
        {
            MatrixBigRational output = new MatrixBigRational(this.Columns, this.Rows);
            for (int i = 0; i < this.Rows; i++)
                for (int j = 0; j < this.Columns; j++)
                    output[j, i] = this[i, j];
            return output;
        }

        private bool GreaterZero()
        {
            for (int i = 0; i < this.Rows; i++)
                for (int j = 0; j < this.Columns; j++)
                    if (this._matrix[i, j] <= 0)
                        return false;
            return true;
        }
    }
}

//        public string ToMathematicaString()
//        {
//            StringBuilder output = new StringBuilder("{");
//            for (int i = 0; i < Rows; i++)
//            {
//                output.Append("{");
//                for (int j = 0; j < Columns; j++)
//                    output.Append(_matrix[i, j]).Append(",");
//                output.Length--;
//                output.Append("},");
//            }
//            output.Length--;
//            output.Append("}");
//            return output.ToString();
//        }

//        public string ToMapleString()
//        {
//            StringBuilder output = new StringBuilder("[");
//            for (int i = 0; i < this.Rows; i++)
//            {
//                output.Append("[");
//                for (int j = 0; j < this.Columns; j++)
//                    output.Append(this[i, j]).Append(",");
//                output.Length--;
//                output.Append("],");
//            }
//            output.Length--;
//            output.Append("]");
//            return output.ToString();
//        }

//        public MatrixBigRational(int dimension, int[] circularLineInteger)
//            //circularLineInteger gibt die Kurzform n/2 oder n/2 + 1 an!!!
//        {
//            BigRational[] matrixArray = new BigRational[dimension];
//            matrixArray[0] = BigRational.Zero;
//            BigRational tempBR;
//
//            if (Functions.IsOdd(dimension))
//            {
//                for (int i = 0; i < circularLineInteger.Length; i++)
//                {
//                    tempBR = new BigRational(new BigInteger(circularLineInteger[i]));
//                    matrixArray[i + 1] = tempBR;
//                    matrixArray[dimension - 1 - i] = tempBR;
//                }
//            }
//            else
//            {
//                for (int i = 0; i < circularLineInteger.Length - 1; i++)
//                {
//                    tempBR = new BigRational(new BigInteger(circularLineInteger[i]));
//                    matrixArray[i + 1] = tempBR;
//                    matrixArray[dimension - 1 - i] = tempBR;
//                }
//                matrixArray[dimension / 2] = new BigRational(new BigInteger(circularLineInteger[circularLineInteger.Length - 1]));
//            }
//
//            Rows =dimension;
//            Columns = dimension;
//            _matrix = new BigRational[Rows, Columns];
//            for (int i = 0; i < Rows; i++)
//            {
//                for (int j = 0; j < Columns; j++)
//                    _matrix[i, j] = matrixArray[j];
//                Functions.ShiftRight(matrixArray);
//            }
//
//            //_rows = circularLineInteger.Length * 2 + 1;
//            //_cols = _rows;
//            //BigRational[] matrixArray = new BigRational[_rows];
//            //matrixArray[0] = BigRational.Zero;
//            //for (int i = 0; i < circularLineInteger.Length; i++)
//            //{
//            //    BigRational tempBR = new BigRational(new BigInteger(circularLineInteger[i]));
//            //    matrixArray[i + 1] = tempBR;
//            //    matrixArray[_rows - 1 - i] = tempBR;
//            //}
//
//            //_matrix = new BigRational[_rows, _cols];
//            //for (int i = 0; i < _rows; i++)
//            //{
//            //    for (int j = 0; j < _cols; j++)
//            //        _matrix[i, j] = matrixArray[j];
//            //    Functions.ShiftRight(matrixArray);
//            //}
//        }
//        public MatrixBigRational(string[,] matrixString)
//        {
//            Rows = matrixString.GetLength(0);
//            Columns = matrixString.GetLength(1);
//            _matrix = new BigRational[Rows, Columns];
//            for (int i = 0; i < Rows; i++)
//                for (int j = 0; j < Columns; j++)
//                    this[i, j] = new BigRational(matrixString[i, j]);
//        }

//        public MatrixBigRational(int[,] matrixInteger)
//        {
//            Rows = matrixInteger.GetLength(0);
//            Columns = matrixInteger.GetLength(1);
//            _matrix = new BigRational[Rows, Columns];
//            for (int i = 0; i < Rows; i++)
//                for (int j = 0; j < Columns; j++)
//                    this[i, j] = new BigRational(new BigInteger(matrixInteger[i, j]));
//        }
//        public MatrixBigRational Clone()
//        {
//            MatrixBigRational mbr = new MatrixBigRational(this.Rows, this.Columns);
//            for (int i = 0; i < mbr.Rows; i++)
//                for (int j = 0; j < mbr.Columns; j++)
//                    mbr[i, j] = this[i, j];
//            return mbr;
//
//        }

//        public string ToMinimalString()
//        {
//            StringBuilder output = new StringBuilder("");
//            for (int i = 0; i < Rows; i++)
//            {
//                output.Append("");
//                for (int j = 0; j < Columns; j++)
//                    output.Append(_matrix[i, j]).Append(",");
//                output.Length--;
//                output.Append(" ");
//            }
//            output.Length--;
//            //output.Append("}");
//            output.Replace("/1,", ",").Replace("/1 ", " ");
//            string outstring = output.ToString();
//            if (outstring.EndsWith("/1"))
//            {
//                output.Length--;
//                output.Length--;
//            }
//            return output.ToString();
//        }

//public static MatrixBigRational MatrixFromCircularDatabaseString(string m)
//{
//    m = m.Replace(" ", "");
//    string[] split = m.Split(new string[] { "," }, StringSplitOptions.RemoveEmptyEntries);
//    int[] intList = new int[split.Length];
//    for (int i = 0; i < split.Length; i++)
//        intList[i] = Convert.ToInt32(split[i]);
//    return new MatrixBigRational(intList);
//}

//public static MatrixBigRational MatrixFromCircularLineIntegerWerner(int[] circularLineInteger)
//{
//    _rows = circularLineInteger.Length * 2 + 1;
//    _cols = _rows;
//    BigRational[] matrixArray = new BigRational[_rows];
//    matrixArray[0] = BigRational.Zero;
//    for (int i = 0; i < circularLineInteger.Length; i++)
//    {
//        BigRational tempBR = new BigRational(new BigInteger(circularLineInteger[i]));
//        matrixArray[i + 1] = tempBR;
//        matrixArray[_rows - 1 - i] = tempBR;
//    }

//    _matrix = new BigRational[_rows, _cols];
//    for (int i = 0; i < _rows; i++)
//    {
//        for (int j = 0; j < _cols; j++)
//            _matrix[i, j] = matrixArray[j];
//        Functions.ShiftRight(matrixArray);
//    }
//}

//        public static MatrixBigRational MatrixFromMathematicaString(string m)
//        {
//            m=m.Replace(" ", "");
//            string[] split = m.Split(new string[] { "},{" }, StringSplitOptions.RemoveEmptyEntries);
//            MatrixBigRational mbr= new MatrixBigRational(split.Length,split[0].Split(new string[] { "," }, StringSplitOptions.RemoveEmptyEntries).Length);
//            for (int i = 0; i < split.Length; i++)
//            {
//                split[i] = split[i].Replace("{", "");
//                split[i] = split[i].Replace("}", "");
//                string[] split2 = split[i].Split(new string[] { "," }, StringSplitOptions.RemoveEmptyEntries);
//                for (int j = 0; j < split2.Length; j++)
//                    mbr[i, j] = new BigRational(split2[j]);
//            }
//            return mbr;
//        }

//        public SquareMatrix ToSquareMatrix()
//        {
//            SquareMatrix output = new SquareMatrix(this._matrix);
//            return output;
//        }

//        public MatrixBigRational GetColAsColVector(int k)
//        {
//            MatrixBigRational m = new MatrixBigRational(Rows, 1);
//            for (int i = 0; i < Rows; i++) 
//                m[i, 0] = _matrix[i, k];
//            return m;
//        }
//
//        public MatrixBigRational GetRowAsColVector(int k)
//        {
//            MatrixBigRational m = new MatrixBigRational(Columns, 1);
//            for (int i = 0; i < Columns; i++)
//                m[i, 0] = _matrix[k, i];
//            return m;
//        }
//
//        public void SetCol(MatrixBigRational v, int k)
//        {
//            for (int i = 0; i < Rows; i++) 
//                _matrix[i, k] = v[i, 0];
//        }
//

//        public MatrixBigRational GetLeadingPrincipalSubmatrix(int length)
//        {
//            MatrixBigRational B = new MatrixBigRational(length);
//            for (int i = 0; i < length; i++)
//                for (int j = 0; j < length; j++)
//                    B[i, j] = this[i, j];
//            return B;
//        }




//public int SignOfDeterminantMathematica()
//{
//    string input = "Sign[Det[" + this.ToMathematicaString() + "]]";
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    return Convert.ToInt32(ret);
//}

//public int SignOfSmallestEigenvalueMathematica()
//{
//    string input = "First[Sort[Sign[Eigenvalues[" + this.ToMathematicaString() + "]]]]";
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    return Convert.ToInt32(ret);
//}

//        public bool IsPositiveDefiniteMathematica()
//        {
//            string input = "PositiveDefiniteMatrixQ[" + this.ToMathematicaString() + "]";
//            string ret = MainClass.MO.EvaluateToInputForm(input);
//            if (ret == "False")
//                return false;
//            else if (ret == "True")
//                return true;
//            else
//                throw new ArgumentException("Mathematica gives an invalid result back", "isPosDef");
//        }

//public bool IsPositiveSemidefiniteMathematica()
//{
//    return (this.SignOfSmallestEigenvalueMathematica() != -1);
//}







//        public MatrixBigRational Vectorize()
//        {
//            MatrixBigRational mbr = new MatrixBigRational(this.Rows * this.Columns, 1);
//            for (int j=0;j<this.Rows;j++)
//                for (int i=0;i<this.Rows;i++)
//                {
//                    mbr[j*this.Rows+i,0]=this[i,j];
//                }
//            return mbr;
//        }

//public int RankMathematica()
//{
//    string mathInput = "Parallelize[MatrixRank[" + this.ToMathematicaString() + "]]";
//    //string ret = Program.MathematicaCall.EnqueueMathTask(new zzzMathematicaInputOutput(mathInput, MathematicaOutputForm.OutputForm));
//    string ret = MainClass.MO.EvaluateToInputForm(mathInput);
//    return Convert.ToInt32(ret);

//}

//public MatrixBigRational (BigRational[] circularLineBigRational)
//{
//    Rows = circularLineBigRational.Length;
//    Columns = Rows;
//    _matrix = new BigRational[Rows, Columns];
//    for (int i = 0; i < Rows; i++)
//    {
//        for (int j = 0; j < Columns; j++)
//            _matrix[i, j] = circularLineBigRational[j];
//        Functions.ShiftRight(circularLineBigRational);
//    }
//}


//        public static MatrixBigRational OneMatrix(int iRows, int iCols)       // Function generates the zero matrix
//        {
//            MatrixBigRational matrix = new MatrixBigRational(iRows, iCols);
//            BigRational brZero = BigRational.One;
//            for (int i = 0; i < iRows; i++)
//                for (int j = 0; j < iCols; j++)
//                    matrix[i, j] = brZero;
//            return matrix;
//        }



//        public static MatrixBigRational ConstantMatrix(int iRows, int iCols, BigRational x)      
//        {
//            MatrixBigRational matrix = new MatrixBigRational(iRows, iCols);
//            for (int i = 0; i < iRows; i++)
//                for (int j = 0; j < iCols; j++)
//                    matrix[i, j] = x;
//            return matrix;
//        }

//        public static MatrixBigRational BlockMatrix(MatrixBigRational A, MatrixBigRational B, MatrixBigRational C, MatrixBigRational D)      
//        {
//            MatrixBigRational matrix = new MatrixBigRational(A.Rows+C.Rows, A.Columns+B.Columns);
//            for (int i = 0; i < matrix.Rows; i++)
//                for (int j = 0; j < matrix.Columns; j++)
//                {
//                    if (i < A.Rows && j < A.Columns)
//                        matrix[i, j] = A[i, j];
//                    else if (i >= A.Rows && j < A.Columns)
//                        matrix[i, j] = C[i - A.Rows, j];
//                    else if (i < A.Rows && j >= A.Columns)
//                        matrix[i, j] = B[i, j - A.Columns];
//                    else if (i >= A.Rows && j >= A.Columns)
//                        matrix[i, j] = D[i - A.Rows, j - A.Columns];
//                }
//            
//            return matrix;
//        }

//        public static MatrixBigRational GetUnitVector(int dimension, int position)
//        {
//            MatrixBigRational mbr = new MatrixBigRational(dimension,1);
//            BigRational zero = BigRational.Zero;
//            for (int i = 0; i < dimension; i++)
//                if (i==position-1)
//                    mbr[i,0] = BigRational.One;
//                else
//                    mbr[i,0] = zero;
//            return mbr;
//        }

//        public static MatrixBigRational MakeMatrixFromListOfVectors(List<MatrixBigRational> vectors)
//        {
//            MatrixBigRational mbr = new MatrixBigRational(vectors.Count, vectors[0].Rows);
//            for (int i = 0; i < mbr.Rows; i++)
//                for (int j = 0; j < mbr.Columns; j++)
//                    mbr[i, j] = vectors[i][j, 0];
//            return mbr;
//        }

//        public static MatrixBigRational CombineTwoEssMatrices (MatrixBigRational A, MatrixBigRational B)
//        {
//            MatrixBigRational mbr = new MatrixBigRational(A.Rows * B.Rows, A.Columns + B.Columns);
//            for (int i_A = 0; i_A < A.Rows; i_A++)
//                for (int i_B = 0; i_B < B.Rows; i_B++)
//                {
//                    for (int j_A = 0; j_A < A.Columns; j_A++)
//                        mbr[A.Rows * i_A + i_B, j_A] = A[i_A, j_A];
//                    for (int j_B = 0; j_B < B.Columns; j_B++)
//                        mbr[A.Rows * i_A + i_B, A.Columns + j_B] = A[i_B, j_B];
//                }
//            return mbr;
//        }

//        public static MatrixBigRational MakeQqtMatrix(MatrixBigRational A)
//        {
//            MatrixBigRational mbr = new MatrixBigRational(A.Columns * A.Columns, A.Rows);
//            for (int j = 0; j < mbr.Columns; j++)
//            {
//                MatrixBigRational x = A.GetRowAsColVector(j);
//                MatrixBigRational qqtVector = (x * (x.Transpose())).Vectorize();
//                x = null;
//                for (int i = 0; i < mbr.Rows; i++)
//                    mbr[i, j] = qqtVector[i, 0];
//
//                qqtVector = null;
//            }
//            return mbr;
//        }

//public static MatrixBigRational MakeQqtMatrix(MatrixBigRational A)
//{
//    MatrixBigRational mbr = new MatrixBigRational(A.Columns * (A.Columns+1)/2, A.Rows);
//    for (int m = 0; m < mbr.Columns; m++)
//        for (int j = 0; j < A.Columns; j++)
//            for (int i = 0; i <=j; i++)
//            { 
//                mbr[j*(j+1)/2+,m]
//            }

//        MatrixBigRational x = A.GetRowAsColVector(j);
//        MatrixBigRational qqtVector = (x * (x.Transpose())).Vectorize();
//        x = null;
//        for (int i = 0; i < mbr.Rows; i++)
//            mbr[i, j] = qqtVector[i, 0];

//        qqtVector = null;
//    }
//    return mbr;
//}
//        public static MatrixBigRational HadamardProduct (MatrixBigRational m1, MatrixBigRational m2)
//        {
//            MatrixBigRational result = new MatrixBigRational(m1.Rows, m1.Columns);
//            for (int i = 0; i < result.Rows; i++)
//                for (int j = 0; j < result.Columns; j++)
//                    result[i, j] = m1[i, j] * m2[i, j];
//            return result;
//        }