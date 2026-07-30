using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace EFR
{
    class EssFinder
    {
        public MatrixBigRational SpielMatrix { get; set; }
        public HashSet<VectorInformation> VectorSet { get; set; }
        public int AnzahlEss { get; set; }
        public string Name { get; set; }
        public long MatrixId { get; set; }
        public int Dimension { get; set; }
        public string Log { get { return _logString.ToString();} }

        private StringBuilder _logString= null;
        private FindEqBigRational _findEquilibriumBR;
        private FindEqDouble _findEqDouble;
        private GrundEss _grundEss;
        private List<int>[] _Jlatein;


        public EssFinder(MatrixBigRational matrix)
        {
            SpielMatrix = matrix;
            Dimension =  matrix.Rows;
            _logString = new StringBuilder("");
            _findEquilibriumBR = new FindEqBigRational(SpielMatrix);
            _findEqDouble = new FindEqDouble(SpielMatrix.ToDoubleArray());
            VectorSet = new HashSet<VectorInformation>();
            if (MainClass.WithLog)
                _logString = new StringBuilder();
        }

        //public string LogToString()
        //{
        //    return _logString.ToString();
        //}

        public void FindESS()
        {
            if (MainClass.WithLog)
            {
                _logString.Append("n=").AppendLine(Dimension.ToString());
                _logString.AppendLine("Spielmatrix:\r\n" + SpielMatrix.ToStringShort());
            }  
            _Jlatein = BitVector.GetAllSubsetsWithMoreThanZeroElementsArrayList(Dimension);

            SearchOneSubsetSize(0);

            for (int i = 1; i <= Dimension - 1; i++)
                SearchOneSubsetSize(i);

            //if (!SearchOneSubsetSize(Dimension - 1))
            //    for (int i = 1; i < Dimension-1; i++)
            //        SearchOneSubsetSize(i);
        }
        

        private bool SearchOneSubsetSize(int JStrichMincount) 
        {
            bool isEss;
            bool hasEss = false;

            int supportSize = JStrichMincount + 1;
            if (MainClass.WithLog)
            {
                _logString.AppendLine("***************************************************************************************************************************************************************************");
                _logString.Append("SupportSize: ").Append(supportSize).AppendLine();
                _logString.Append("Anzahl Jmin:").AppendLine(_Jlatein[JStrichMincount].Count.ToString());
            }

            //find equilibria and stability from jmin
            foreach (int elementOfJmin in _Jlatein[JStrichMincount])
            {
                if (MainClass.WithLog)
                    _logString.Append("FINDEQ aufgerufen für: ").AppendLine(BitVector.ToListOfInt(elementOfJmin).ToStringExtended());

                bool ResultsFindeq = false;
                if (_findEqDouble.IsFindEq(elementOfJmin, supportSize))
                    ResultsFindeq = _findEquilibriumBR.FindEquilibrium(elementOfJmin);//, supportSize, out extSupport, out extSupportDimension, out payoff);

                if (ResultsFindeq)
                {
                    //überprüfe checkstab
                    isEss = false;
                    if (CHECKSTAB(_findEquilibriumBR.ResultVector, elementOfJmin, _findEquilibriumBR.ExtendedSupport, _findEquilibriumBR.ExtendedSupportSize))
                    {
                        isEss = true;
                        hasEss = true;
                        AnzahlEss += 1;
                    }
                    if (MainClass.WithLog)
                        _logString.Append("CHECKSTAB: ").AppendLine(_grundEss.ToString());

                    VectorInformation item = new VectorInformation();
                    item.Support = BitVector.ToListString(elementOfJmin);
                    item.SupportSize = supportSize;
                    item.ExtendedSupport = BitVector.ToListString(_findEquilibriumBR.ExtendedSupport);
                    item.ExtendedSupportSize = _findEquilibriumBR.ExtendedSupportSize;
                    item.IsEss = isEss;
                    item.Vektor = _findEquilibriumBR.ResultVector.ToRString();
                    item.GrundEss = (byte)_grundEss;
                    item.Payoff = (decimal)_findEquilibriumBR.Payoff;
                    VectorSet.Add(item);

                    //entferne alle sets aus jstrich für diesen support
                    for (int i = JStrichMincount + 1; i < Dimension; i++)
                        _Jlatein[i].RemoveAll(x => (elementOfJmin & x) == elementOfJmin);
                }
            }
            return hasEss;
        }

        private bool CHECKSTAB(MatrixBigRational p, int IofP, int JofP, int JofPCount)
        {
            int m = (IofP & (~(IofP - 1))); //get lowest set bit
            int indexOfMreal=BitVector.GetPostitionOfIthSetFlag(m,1)-1;
            int JStrichofP = JofP & (~m); //ext support without m

            if (MainClass.WithLog)
            {
                _logString.AppendLine("---------------------------------------------------------------------------------------------------------------------");
                _logString.AppendLine("CHECKSTAB aufgerufen für: " + p.ToMathematicaString());
                _logString.AppendLine("Iofp: " + BitVector.ToListOfInt( IofP).ToStringExtended());
                _logString.AppendLine("Jofp: " + BitVector.ToListOfInt(JofP).ToStringExtended());
                _logString.AppendLine("m: " + BitVector.ToListOfInt(m).ToStringExtended());
                _logString.AppendLine("JStrichofp: " + BitVector.ToListOfInt(JStrichofP).ToStringExtended());
            }
            int JStrichofPCount = JofPCount - 1;
            if (JStrichofPCount == 0) 
            {
                _grundEss = GrundEss.true_pureEss;
                return true;
            }

            SquareMatrix B = new SquareMatrix(JStrichofPCount);

            int row = 0;
            int column = 0;
            for (int i = 0;i<Dimension;i++)
                if ((JStrichofP & (1 << i)) != 0)
                {
                    column = 0;
                    for (int j = 0; j < Dimension; j++)
                        if ((JStrichofP & (1 << j)) != 0)
                        {
                            BigRational x = SpielMatrix[indexOfMreal, j] + SpielMatrix[j, indexOfMreal] +
                                SpielMatrix[i, indexOfMreal] + SpielMatrix[indexOfMreal, i] -
                                SpielMatrix[i, j] - SpielMatrix[j, i] - 2 * SpielMatrix[indexOfMreal, indexOfMreal];
                            B.SetElement(row,column,i+1,j+1, x);
                            column += 1;
                        }
                    row += 1;
                }
            MatrixBigRational BBR = B.ToMatrixBigRational();
            if (BBR.IsPosDefDouble())
            {
                _grundEss = GrundEss.true_posdef_double;
                return true;
            }
            int K = (JofP & (~IofP)); //JofP except IofP
            int Ksize = BitVector.SupportSize(K);
            if (Ksize==0)
            {
                if (BBR.IsPositiveDefiniteMathematica())
                {
                    _grundEss = GrundEss.true_copos_BR;
                    return true;
                }
                else
                {
                    _grundEss = GrundEss.false_copos_BR;
                    return false;
                }
            }
            if (IsCoposWithRespectToK(B, BitVector.ToListOfInt(JStrichofP), BitVector.ToListOfInt(K)))
                return true;
            else
                return false;
        }

        private bool IsCoposWithRespectToK(SquareMatrix B, List<int> J, List<int> K)
        {
            //überprüfe zuerst, ob "einfache" lösungen bestehen
            if (MainClass.WithLog)
            {
                _logString.AppendLine("Copos with respect to aufgerufen:");
                _logString.AppendLine("Menge J: " + J.ToStringExtended());
                _logString.AppendLine("Menge K: " + K.ToStringExtended());
                _logString.AppendLine("Matrix B:").AppendLine(B.ToMatrixBigRational().ToMathematicaString());
            }

            ////bomze thr 3.4 a
            if (K.Count == 2 && J.Count == 2)
            {
                if ((B.GetValue(J[0], J[0]) > 0 && B.GetValue(J[1], J[1]) > 0) &&
                            ((B.GetValue(J[0], J[0]) * B.GetValue(J[1], J[1]) - B.GetValue(J[0], J[1]) * B.GetValue(J[1], J[0])) > 0 || B.GetValue(J[0], J[1]) >= 0))
                {
                    _grundEss = GrundEss.true_copos_wrt_Kequals2andJequals2;
                    return true;
                }  
                else
                {
                    _grundEss = GrundEss.false_copos_wrt_Kequals2andJequals2;
                    return false;
                }
            }

            //keine einfachen lösungen, rekursive copositivity erforderlich!
            List<int> JohneK = new List<int>(J.Except(K));
            List<List<int>> listK = new List<List<int>>();
            List<SquareMatrix> listB = new List<SquareMatrix>();
            JohneK = J.Except(K).ToList();
            listK.Add(J);
            listB.Add(B);
            int n;

            for (int v = 1; v <= JohneK.Count; v++)
            {
                BigRational diagElem = listB[v - 1].GetValue(JohneK[v - 1], JohneK[v - 1]);
                if (diagElem <= BigRational.Zero)
                {
                    _grundEss = GrundEss.false_copos_wrt_diag_leq_0;
                    return false;
                }
                listK.Add(new List<int>(listK[v - 1].Except(Functions.ListofElement(JohneK[v - 1]))));
                n = listK[v].Count;
                listB.Add(new SquareMatrix(n));
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                        listB[v].SetElement(j, k, listK[v][j], listK[v][k], listB[v - 1].GetValue(JohneK[v - 1], JohneK[v - 1]) * listB[v - 1].GetValue(listK[v][j], listK[v][k]) - listB[v - 1].GetValue(JohneK[v - 1], listK[v][j]) * listB[v - 1].GetValue(JohneK[v - 1], listK[v][k]));
            }

            if (K.Count == 0)
            {
                _grundEss = GrundEss.true_copos_wrt_Kequals0_posdef;
                return true;
            }
            //if (listB[JohneK.Count].IsCopositive(K))
            if (listB[JohneK.Count].ToMatrixBigRational().IsCopositive())
            {
                _grundEss = GrundEss.true_copos;
                return true;
            } 
            else
            {
                _grundEss = GrundEss.false_copos;
                return false;
            }
        }
    }
}