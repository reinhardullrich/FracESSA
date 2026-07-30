using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;



namespace EFR
{
    class SquareMatrix
    {
        private Element[,] A;
        public int n { get; private set; }

        struct Element
        {
            public int Row { get; set; }
            public int Column { get; set; }
            public BigRational Value { get; set; }
            public Element(int i, int j, BigRational x)
                : this()
            {
                Row = i;
                Column = j;
                Value = x;
            }
        };

        public SquareMatrix(int n)
        {
            A = new Element[n, n];
            this.n = n;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = new Element(-999, -999, BigRational.MinusOne);
                }
            }
        }

        public SquareMatrix(BigRational[,] varray)
        {
            if (varray.GetLength(0) != varray.GetLength(1))
                throw new ArgumentException("The matrix is not a square matrix", "varray");
            n = varray.GetLength(0);
            A = new Element[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    A[i, j] = new Element(i + 1, j + 1, varray[i, j]);
        }

        public void SetElement(int i, int j, int row, int column, BigRational value)
        {
            A[i, j].Row = row;
            A[i, j].Column = column;
            A[i, j].Value = value;

        }

        public BigRational GetValue(int row, int column)
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (A[i, j].Row == row && A[i, j].Column == column)
                        return A[i, j].Value;
            throw new ArgumentException("The Element(" + row + "," + column + ") is not in the boundaries of the matrix", "row/column");
        }


        public MatrixBigRational ToMatrixBigRational()
        {
            MatrixBigRational mbr = new MatrixBigRational(n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    mbr[i, j] = A[i, j].Value;

            return mbr;
        }





        //public BigRational GetValueByRealIndex(int realRow, int realColumn)
        //{
        //    return A[realRow, realColumn].Value;
        //}



        //public SquareMatrix(int[,] varray)
        //{
        //    if (varray.GetLength(0) != varray.GetLength(1))
        //        throw new ArgumentException("The matrix is not a square matrix", "varray");
        //    n = varray.GetLength(0);
        //    A = new Element[n, n];
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //        for (int j = 0; j < n; j++)
        //            A[shiftDim, j] = new Element(shiftDim + 1, j + 1, new BigRational(Convert.ToDecimal(varray[shiftDim, j])));
        //}

        //public SquareMatrix(string[,] varray)
        //{
        //    if (varray.GetLength(0) != varray.GetLength(1))
        //        throw new ArgumentException("The matrix is not a square matrix", "varray");
        //    n =  varray.GetLength(0);
        //    A = new Element[n, n];
        //    for (int i = 0; i < n; i++)
        //        for (int j = 0; j < n; j++)
        //            try
        //            {
        //                A[i, j] = new Element(i + 1, j + 1, new BigRational(varray[i, j]));
        //            }
        //            catch 
        //            {
        //                throw new ArgumentException("The Element(" + i + "," + j + ") is not valid fraction", "varray");
        //            }
        //}


        //public SquareMatrix(int[,] intArray)
        //{
        //    if (intArray.GetLength(0) != intArray.GetLength(1))
        //        throw new ArgumentException("The matrix is not a square matrix", "varray");
        //    n = intArray.GetLength(0);
        //    A = new Element[n, n];
        //    for (int i = 0; i < n; i++)
        //        for (int j = 0; j < n; j++)
        //            A[i, j] = new Element(i + 1, j + 1, new BigRational(new BigInteger(intArray[i, j])));
        //}

        //public SquareMatrix(int[] intArray)
        //{
        //    this.n = intArray.Length*2+1;
        //    this.A = new Element[n, n];

        //    BigRational[] circularLine = new BigRational[n];

        //    circularLine[0] = BigRational.Zero;
        //    for (int i = 0; i < intArray.Length; i++)
        //    {
        //        BigRational tempBR = new BigRational(new BigInteger(intArray[i]));
        //        circularLine[i + 1] = tempBR;
        //        circularLine[n-1 - i] = tempBR;
        //    }

        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < n; j++)
        //            A[i, j] = new Element(i + 1, j + 1, circularLine[j]);
        //        Functions.ShiftRight(circularLine);
        //    }

        //}

        //public SquareMatrix(SquareMatrix eMatrix, Vector column, Vector row, BigRational number)
        //{
        //    this.n = eMatrix.n + 1;
        //    A = new Element[n, n];
        //    for (int shiftDim = 0; shiftDim < n - 1; shiftDim++)
        //        for (int j = 0; j < n - 1; j++)
        //            A[shiftDim, j] =  new Element(shiftDim + 1, j + 1, eMatrix.GetValueByRealIndex(shiftDim, j));
        //    for (int shiftDim=0;shiftDim<n-1;shiftDim++)
        //        A[shiftDim,n-1]= new Element(shiftDim+1,n,column.GetValue(shiftDim));
        //    for (int j = 0; j < n - 1; j++)
        //        A[n-1, j] = new Element(n, j+1, row.GetValue(j));
        //    A[n - 1, n - 1] = new Element(n, n, number);
        //}

        //public SquareMatrix (BigRational[] circularLine)
        //{
        //    this.n=circularLine.Length;
        //    A= new Element[n,n];
        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < n; j++)
        //            A[i,j]= new Element(i+1,j+1,circularLine[j]);
        //        Functions.ShiftRight(circularLine);
        //    }
        //}

        //public SquareMatrix(string[] circularLine)
        //{
        //    this.n = circularLine.Length;
        //    BigRational[] brLine = new BigRational[n];
        //    for (int i=0;i<n;i++)
        //        brLine[i]=new BigRational(circularLine[i]);
        //    A = new Element[n, n];
        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < n; j++)
        //            A[i, j] = new Element(i + 1, j + 1, brLine[j]);
        //        Functions.ShiftRight(brLine);
        //    }

        //}



        //public SquareMatrix Transpose()
        //{
        //    SquareMatrix outMatrix = new SquareMatrix(this.n);
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //        for (int j = 0; j < n; j++)
        //            outMatrix.SetElement(shiftDim,j,shiftDim+1,j+1, this.GetValueByRealIndex(j, shiftDim));
        //    return outMatrix;
        //}


        //public Vector GetColumn(int column)
        //{
        //    Vector returnVector = new Vector(n);
        //    for (int shiftDim=0;shiftDim<n;shiftDim++)
        //        returnVector.SetValue(shiftDim,A[shiftDim, column].Value);
        //    return returnVector;
        //}

        //public Vector GetRow(int row)
        //{
        //    Vector returnVector = new Vector(n);
        //    for (int shiftDim=0;shiftDim<n;shiftDim++)
        //        returnVector.SetValue(shiftDim,A[row, shiftDim].Value);
        //    return returnVector;
        //}

        //public override string ToString()
        //{
        //    string str = "";
        //    for (int i = 0; i < A.GetLength(0); i++)
        //    {
        //        for (int j = 0; j < A.GetLength(1); j++)
        //            str += "a(" + A[i, j].Row + "," + A[i, j].Column + ")= " + A[i, j].Value + "  ";
        //        str += "\r\n";
        //    }
        //    return str.Left(str.Length - 4);
        //}

        //public  string ToStringShort()
        //{
        //    string str = "";
        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int j = 0; j < n; j++)
        //            str +=  A[i, j].Value + "\t";
        //        str += "\r\n";
        //    }
        //    return str.Left(str.Length - 4);
        //}

        //public BigRational GetQuadraticForm(Vector p)
        //{
        //    if (p.Dimension!=n)
        //        throw new ArgumentException("Not the same dimension","p");
        //    Vector Ap = this.MultiplyWithVector(p);
        //    BigRational outPut = Ap.Dot(p);
        //    return outPut;
        //}

        //public SquareMatrix GetSubmatrix(List<int> indexer)
        //{
        //    SquareMatrix B = new SquareMatrix(indexer.Count);
        //    for (int shiftDim = 0; shiftDim < indexer.Count; shiftDim++)
        //        for (int j = 0; j < indexer.Count; j++)
        //            B.SetElement(shiftDim, j, indexer[shiftDim], indexer[j], this.GetValue(indexer[shiftDim], indexer[j]));
        //    return B;
        //}

        //public Vector MultiplyWithVector(Vector p)
        //{
        //    if (p.Dimension != n)
        //        throw new ArgumentException("The dimensions of the matrix and the vector are not equal", "vector");
        //    Vector returnVector = new Vector(n);
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //    {
        //        BigRational rowSum=BigRational.Zero;
        //        for (int j = 0; j < n; j++)
        //            rowSum += A[shiftDim, j].Value * p.GetValue(j);
        //        returnVector.SetValue(shiftDim, rowSum);
        //    }
        //    return returnVector;
        //}

        //private ILArray<double> ToLNumericsMatrix()
        //{
        //    ILArray<double> matrixLNum = ILMath.zeros(n, n);
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //        for (int j = 0; j < n; j++)
        //            matrixLNum[shiftDim, j] = (double)this.A[shiftDim, j].Value;
        //    return matrixLNum;
        //}

        //public string ToMathMatrix()
        //{
        //    StringBuilder output = new StringBuilder("{");
        //    for (int i = 0; i < n; i++)
        //    {
        //        output.Append("{");
        //        for (int j = 0; j < n; j++)
        //            output.Append(this.A[i, j].Value).Append(",");
        //        output.Length--;
        //        output.Append("},");
        //    }
        //    output.Length--;
        //    output.Append("}");
        //    return output.ToString();
        //}

        //public string ToMapleMatrix()
        //{
        //    StringBuilder output = new StringBuilder("[");
        //    for (int i = 0; i < n; i++)
        //    {
        //        output.Append("[");
        //        for (int j = 0; j < n; j++)
        //            output.Append(this.A[i, j].Value).Append(",");
        //        output.Length--;
        //        output.Append("],");
        //    }
        //    output.Length--;
        //    output.Append("]");
        //    return output.ToString();
        //}

        //public decimal[,] ToDecimalArray()
        //{
        //    decimal[,] output = new decimal[n, n];
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //        for (int j = 0; j < n; j++)
        //            output[shiftDim, j] = ((decimal)this.A[shiftDim, j].Value);
        //    return output;
        //}

        //public double[,] ToDoubleArray()
        //{
        //    double[,] output = new double[n, n];
        //    for (int i = 0; i < n; i++)
        //        for (int j = 0; j < n; j++)
        //            output[i, j] = ((double)this.A[i, j].Value);
        //    return output;
        //}

        //public string[,] ToStringArray()
        //{
        //    string[,] ret = new string[n, n];
        //    for (int shiftDim = 0; shiftDim < n; shiftDim++)
        //        for (int j = 0; j < n; j++)
        //            ret[shiftDim, j] = this.A[shiftDim, j].Value.ToString()+",";
        //    return ret;
        //}

        //public BigRational[,] ToBigRationalArray()
        //{
        //    BigRational [,] output = new BigRational [n, n];
        //    for (int i = 0; i < n; i++)
        //        for (int j = 0; j < n; j++)
        //            output[i, j] = this.A[i, j].Value;
        //    return output;
        //}

        //public bool IsPosDef()
        //{
        //    string output="";
        //    output = Program.MathOrg.EvaluateToInputForm("PositiveDefiniteMatrixQ[" + this.ToMathMatrix() + "]");
        //    if (output == "False")
        //        return false;
        //    else if (output == "True")
        //        return true;
        //    else
        //        throw new ArgumentException("Mathematica gives an invalid result back", "IsPosDef");


        //ILArray<double> matrixLNum = this.ToLNumericsMatrix();
        //try
        //{
        //    ILMath.chol(matrixLNum, true);
        //}
        //catch (ILNumerics.Exceptions.ILArgumentException ex) //throws exception if not pos.def.
        //{
        //    return false; 
        //}
        //return true;

        //}

        //        public bool IsCopositive(List<int> K)
        //        {
        //#if DEBUG
        //            //System.Diagnostics.Debugger.Break();
        //#endif

        //            bool coPos = true;
        //            LinkedList<List<int>> PowerSet = Functions.PowerSet(K, 0);
        //            foreach (List<int> Set in PowerSet)
        //            {
        //                SquareMatrix C = this.GetSubmatrix(Set);
        //                ILArray<double> C_il = C.ToLNumericsMatrix();
        //                ILArray<double> EV = ILMath.ones(C_il.Count);
        //                ILArray<double> EW = ILMath.eigSymm(C_il, EV);
        //                Vector ev;
        //                for (int shiftDim = 0; shiftDim < EV.Length; shiftDim++)
        //                {
        //                    ev = new Vector(EV[ILMath.full, shiftDim]);
        //                    if (ev.GetMinValue() > BigRational.Zero)
        //                        if (new BigRational((double)EW[shiftDim, shiftDim]) <= 0)
        //                        {
        //                            coPos = false;
        //                            break;
        //                        }
        //                } 
        //            }
        //            return coPos;
        //        }

        //public static SquareMatrix operator +(SquareMatrix left, SquareMatrix right)
        //{
        //    if (left.n != right.n)
        //        throw new ArgumentException("Matrices don't have the same dimension", "left/right");

        //    SquareMatrix outMatrix = new SquareMatrix(left.n);
        //    for (int i = 0; i < left.n; i++)
        //        for (int j = 0; j < left.n; j++)
        //            outMatrix.SetElement(i, j, i + 1, j + 1, left.GetValueByRealIndex(j, i) + right.GetValueByRealIndex(j, i));
        //    return outMatrix;
        //}

        //public  SquareMatrix GetLeadingPrincipalSubmatrix( int length)
        //{
        //    SquareMatrix B = new SquareMatrix(length);
        //    for (int i = 0; i < length; i++)
        //        for (int j = 0; j < length; j++)
        //            B.SetElement(i, j, i + 1, j + 1, this.GetValueByRealIndex(i, j));
        //    return B;
        //}


        //public static SquareMatrix BlockTwoMatrixWithConstants( SquareMatrix leftBlock, SquareMatrix rightBlock, BigRational leftValue, BigRational rightValue)
        //{
        //    int n = leftBlock.n + rightBlock.n;
        //    SquareMatrix A = new SquareMatrix(n);
        //    for (int i =0;i<n; i++)
        //        for (int j=0;j<n;j++)
        //            if (i<leftBlock.n)
        //                if (j<leftBlock.n)
        //                    A.SetElement(i,j,i+1,j+1,leftBlock.GetValueByRealIndex(i,j));
        //                else
        //                    A.SetElement(i,j,i+1,j+1,rightValue);
        //            else
        //                if (j<leftBlock.n)
        //                    A.SetElement(i,j,i+1,j+1,leftValue);
        //                else
        //                    A.SetElement(i,j,i+1,j+1,rightBlock.GetValueByRealIndex(i-leftBlock.n,j-leftBlock.n));
        //    return A;
        //}

        //public static SquareMatrix IdentityMatrix(int n)
        //{
        //    SquareMatrix A = new SquareMatrix(n);
        //    BigRational one = BigRational.One;
        //    BigRational zero=BigRational.Zero;
        //    for (int i =0;i<n; i++)
        //        for (int j=0;j<n;j++)
        //            if (i==j)
        //                A.SetElement(i,j,i+1,j+1,one);
        //            else
        //                A.SetElement(i,j,i+1,j+1,zero);
        //    return A;
        //}




    }
}

