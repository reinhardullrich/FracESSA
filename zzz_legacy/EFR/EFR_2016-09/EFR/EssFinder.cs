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
        public HashSet<VectorInformation> VectorSet { get; private set; }
        public int NumberOfEss { get; private set; }
		public int NumberOfCandidates { get; private set; }
        public string Name { get; set; }
        public long MatrixId { get; set; }
        public int Dimension { get; private set; }
        public string Log { get { return _log.ToString();} }
		public bool IsIndeterminate { get; private set;}

        private StringBuilder _log;
		private FindNashRational _findNashRational;
		private FindNashDouble _findNashDouble;
        private List<int>[] _Jlatein;

        public EssFinder(MatrixBigRational matrix)
        {
            SpielMatrix = matrix;
            Dimension =  matrix.Rows;
			IsIndeterminate = false;
			VectorSet = new HashSet<VectorInformation>();
            _log = new StringBuilder("");
            _findNashRational = new FindNashRational(SpielMatrix);
            _findNashDouble = new FindNashDouble(SpielMatrix.ToDoubleArray());
            VectorSet = new HashSet<VectorInformation>();
            if (EFR.WithLog)
				_log = new StringBuilder();
        }
			
		public string GetVectorCsv ()
		{
			StringBuilder x = new StringBuilder ("");
			x.AppendLine (VectorInformation.GetHeader ());
			foreach (VectorInformation vector in this.VectorSet)
			{  
				x.AppendLine (vector.GetVectorInformation());            
			}
			x.Length--;
			return x.ToString();
		}

        public void FindESS()
        {
            if (EFR.WithLog) {
                _log.Append("n=").AppendLine(Dimension.ToString());
                _log.AppendLine("Gamematrix:\r\n" + SpielMatrix.ToMatrixRepresenation());
            }  
            _Jlatein = BitVector.GetAllSubsetsWithMoreThanZeroElementsArrayList(Dimension);
            for (int i = 0; i < Dimension; i++)
                SearchOneSubsetSize(i);
			if (IsIndeterminate)
				NumberOfEss = -1;
			if (EFR.WithLog) {
				_log.AppendLine ("Number of isolated NES: " + this.NumberOfCandidates);
				_log.AppendLine ("Number of ESS: " + this.NumberOfEss);
			}
				
        }
        
        private void SearchOneSubsetSize(int JStrichMincount) 
        {
			int checkstabResult;
			ReasonEss reasonEss;

            int supportSize = JStrichMincount + 1;

            if (EFR.WithLog)
            {
                _log.AppendLine("***************************************************************************************************************************************************************************");
                _log.Append("SupportSize: ").Append(supportSize).AppendLine();
                _log.Append("Number of Elements in Jmin:").AppendLine(_Jlatein[JStrichMincount].Count.ToString());
            }

            //find equilibria and stability from jmin
            foreach (int elementOfJmin in _Jlatein[JStrichMincount])
            {
                if (EFR.WithLog)
                    _log.Append("FINDEQ for: ").AppendLine(BitVector.ToListOfInt(elementOfJmin).ToStringExtended());

				if (!_findNashDouble.IsNashDouble (elementOfJmin, supportSize))
					continue;				
				if (!_findNashRational.GetNashRational (elementOfJmin, supportSize))
					continue;

                //run checkstab

				checkstabResult = CHECKSTAB(_findNashRational.ResultVector, elementOfJmin, _findNashRational.ExtendedSupport, _findNashRational.ExtendedSupportSize, out reasonEss);
                
				NumberOfCandidates += 1;

				if (checkstabResult==1) {
                    NumberOfEss += 1;
                }
				else if (checkstabResult==-1) {
					IsIndeterminate = true;
				}
                if (EFR.WithLog)
                    _log.Append("CHECKSTAB: ").AppendLine(reasonEss.ToString());

                VectorInformation item = new VectorInformation();
				item.VectorID = NumberOfCandidates;
                item.Support = BitVector.ToListString(elementOfJmin);
                item.SupportSize = supportSize;
                item.ExtendedSupport = BitVector.ToListString(_findNashRational.ExtendedSupport);
                item.ExtendedSupportSize = _findNashRational.ExtendedSupportSize;
                item.IsEss = checkstabResult;
                item.Vector = _findNashRational.ResultVector.ToString();
				item.ReasonEss = (int)reasonEss;
				item.Payoff = _findNashRational.Payoff;
                item.PayoffDecimal = (decimal)_findNashRational.Payoff;
                VectorSet.Add(item);

                //entferne alle sets aus jstrich für diesen support
                for (int i = JStrichMincount + 1; i < Dimension; i++)
                    _Jlatein[i].RemoveAll(x => (elementOfJmin & x) == elementOfJmin);
            }
        }

		private int CHECKSTAB(MatrixBigRational p, int IofP, int JofP, int JofPCount, out ReasonEss reasonEss)
        {
            int bitsetm = (IofP & (~(IofP - 1))); //get lowest set bit
			int m=BitVector.GetSmallestSetFlag(bitsetm); //get position (starting at 0) for the rightmost set bit
            int JStrichofP = JofP & (~bitsetm); //ext support without m

            if (EFR.WithLog)
            {
                _log.AppendLine("---------------------------------------------------------------------------------------------------------------------");
                _log.AppendLine("CHECKSTAB for: " + p.ToString());
                _log.AppendLine("Iofp: " + BitVector.ToListOfInt( IofP).ToStringExtended());
                _log.AppendLine("Jofp: " + BitVector.ToListOfInt(JofP).ToStringExtended());
                _log.AppendLine("m: " + BitVector.ToListOfInt(bitsetm).ToStringExtended());
                _log.AppendLine("JStrichofp: " + BitVector.ToListOfInt(JStrichofP).ToStringExtended());
            }
            int JStrichofPCount = JofPCount - 1;
            if (JStrichofPCount == 0) 
            {
				reasonEss = ReasonEss.true_pure_ess;
                return 1;
            }

            MatrixBigRational B = new MatrixBigRational(JStrichofPCount);

            int row = 0;
            int column = 0;
            for (int i = 0;i<Dimension;i++)
                if ((JStrichofP & (1 << i)) != 0) {
                    column = 0;
					for (int j = 0; j < i+1; j++)
                        if ((JStrichofP & (1 << j)) != 0) {
							B[row,column] = B[column,row] = SpielMatrix[m, j] + SpielMatrix[j, m] + SpielMatrix[i, m] + SpielMatrix[m, i] -
                                SpielMatrix[i, j] - SpielMatrix[j, i] - 2 * SpielMatrix[m, m];
                            column += 1;
                        }
                    row += 1;
                }

			if (B.IsPosDefDouble()) {
				reasonEss = ReasonEss.true_posdef_double;
                return 1;
            }
				
			if (B.IsPosDef()) {
				reasonEss = ReasonEss.true_posdef_rational;
				return 1;
			}

            int K = (JofP & (~IofP)); //JofP except IofP
            int Ksize = BitVector.SupportSize(K);
			if (Ksize==0 || Ksize==1) {
				reasonEss = ReasonEss.false_not_posdef_and_K_0_1;
	            return 0;
            }
			reasonEss = ReasonEss.indeterminate;
			return -1;
        }
    }
}


//            if (IsCoposWithRespectToK(B, BitVector.ToListOfInt(JStrichofP), BitVector.ToListOfInt(K)))
//                return true;
//            else
//                return false;

//        private bool IsCoposWithRespectToK(SquareMatrix B, List<int> J, List<int> K)
//        {
//            //überprüfe zuerst, ob "einfache" lösungen bestehen
//            if (EFR.WithLog)
//            {
//                _logString.AppendLine("Copos with respect to aufgerufen:");
//                _logString.AppendLine("Menge J: " + J.ToStringExtended());
//                _logString.AppendLine("Menge K: " + K.ToStringExtended());
//                _logString.AppendLine("Matrix B:").AppendLine(B.ToMatrixBigRational().ToMathematicaString());
//            }
//
//            ////bomze thr 3.4 a
//            if (K.Count == 2 && J.Count == 2)
//            {
//                if ((B.GetValue(J[0], J[0]) > 0 && B.GetValue(J[1], J[1]) > 0) &&
//                            ((B.GetValue(J[0], J[0]) * B.GetValue(J[1], J[1]) - B.GetValue(J[0], J[1]) * B.GetValue(J[1], J[0])) > 0 || B.GetValue(J[0], J[1]) >= 0))
//                {
//                    _grundEss = GrundEss.true_copos_wrt_Kequals2andJequals2;
//                    return true;
//                }  
//                else
//                {
//                    _grundEss = GrundEss.false_copos_wrt_Kequals2andJequals2;
//                    return false;
//                }
//            }
//
//            //keine einfachen lösungen, rekursive copositivity erforderlich!
//            List<int> JohneK = new List<int>(J.Except(K));
//            List<List<int>> listK = new List<List<int>>();
//            List<SquareMatrix> listB = new List<SquareMatrix>();
//            JohneK = J.Except(K).ToList();
//            listK.Add(J);
//            listB.Add(B);
//            int n;
//
//            for (int v = 1; v <= JohneK.Count; v++)
//            {
//                BigRational diagElem = listB[v - 1].GetValue(JohneK[v - 1], JohneK[v - 1]);
//                if (diagElem <= BigRational.Zero)
//                {
//                    _grundEss = GrundEss.false_copos_wrt_diag_leq_0;
//                    return false;
//                }
//                listK.Add(new List<int>(listK[v - 1].Except(Functions.ListofElement(JohneK[v - 1]))));
//                n = listK[v].Count;
//                listB.Add(new SquareMatrix(n));
//                for (int j = 0; j < n; j++)
//                    for (int k = 0; k < n; k++)
//                        listB[v].SetElement(j, k, listK[v][j], listK[v][k], listB[v - 1].GetValue(JohneK[v - 1], JohneK[v - 1]) * listB[v - 1].GetValue(listK[v][j], listK[v][k]) - listB[v - 1].GetValue(JohneK[v - 1], listK[v][j]) * listB[v - 1].GetValue(JohneK[v - 1], listK[v][k]));
//            }
//
//            if (K.Count == 0)
//            {
//                _grundEss = GrundEss.true_copos_wrt_Kequals0_posdef;
//                return true;
//            }
//            //if (listB[JohneK.Count].IsCopositive(K))
//            if (listB[JohneK.Count].ToMatrixBigRational().IsCopositive())
//            {
//                _grundEss = GrundEss.true_copos;
//                return true;
//            } 
//            else
//            {
//                _grundEss = GrundEss.false_copos;
//                return false;
//            }
//        }