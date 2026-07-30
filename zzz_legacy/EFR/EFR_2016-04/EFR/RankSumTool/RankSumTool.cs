//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace EssFinderProject.RankSumTool
//{
//    class RankSumTool
//    {

//        private const int MaxDimension=100;
//        private const int MaxIterations = 5;

//        private struct DimensionInfo
//        {
//            public int Dimension = -1;
//            public int d=-1;
//            public int s = -1;
//            public LinkedList<RankInfo> RankLL = new LinkedList<RankInfo>();
//            public int MaxCpRank = -1;

//            public DimensionInfo(int dimension, int d, int s, int maxCpRank)
//            {
//                this.Dimension = dimension;
//                this.d = d;
//                this.s = s;
//                this.MaxCpRank = maxCpRank;
//            }
//        }
//        private struct RankInfo
//        {
//            public int Dimension = -1;
//            public int Rank = -1;
//            public int Rank2=-1;
//            public string Herkunft = "xxx";
//            public int NumberOfBaseMatrices = -1;


//            public RankInfo( int dimension, int rank, int rank2, string herkunft, int numberOfBaseMatrices)
//            {
//                this.Dimension = dimension;
//                this.Rank = rank;
//                this.Rank2 = rank2;
//                this.Herkunft = herkunft;
//                this.NumberOfBaseMatrices = numberOfBaseMatrices;
//            }
//        }

//        private DimensionInfo[,] _DimensionArray = new DimensionInfo[MaxIterations, MaxDimension];

//        public RankSumTool()
//        {
//            Init();

//            for (int i = 1; i <= MaxIterations; i++)
//            {
//                for (int j = 0; j < MaxDimension; j++)
//                {
//                    _RankArray[i, j] = _RankArray[i - 1, j];
//                }
//                for (int j = 0; j < MaxDimension; j++)
//                    for (int k = 0; k < MaxDimension; k++)
//                    {
//                        _RankArray[i, j] = _RankArray[i - 1, j];


//                    }

//            }
                

//        }




//        private void Init()
//        {
//            for (int i = 0; i < MaxDimension;i++ )
//            {
//                _DimensionArray[0, i] = new DimensionInfo(i + 1, d_n(i + 1), s_n(i + 1),i+1);
//                RankInfo rankInfo = new RankInfo( i + 1, i + 1, i + 1, "I_" + (i + 1),1);
//                _DimensionArray[0, i].RankLL.AddLast(rankInfo);
                
//            }

//            _DimensionArray[0, 5] = new RankInfo(6, 6, 8, "U_6", d_n(6), s_n(6));
//            _DimensionArray[0, 6] = new RankInfo(7, 7, 14, "U_7", d_n(7), s_n(7));
//            _DimensionArray[0, 7] = new RankInfo(8, 8, 18, "U_8", d_n(8), s_n(8));
//            _DimensionArray[0, 8] = new RankInfo(9, 9, 26, "U_9", d_n(9), s_n(9));
//            _DimensionArray[0, 10] = new RankInfo(11, 11, 32, "U_11", d_n(11), s_n(11));
//            _DimensionArray[0, 12] = new RankInfo(13, 13, 50, "U_13", d_n(13), s_n(13));
//            _DimensionArray[0, 14] = new RankInfo(15, 15, 360, "U_15", d_n(15), s_n(15));
//        }


//        private int d_n(int dimension)
//        {
//            return (int)Math.Floor((decimal)dimension * (decimal)dimension / 4);
//        }

//        private int s_n(int dimension)
//        {
//            return dimension * (dimension + 1) / 2 - 4;
//        }


//        private RankInfo Calculate(RankInfo x, RankInfo y)
//        {
//            RankInfo ret = new RankInfo();
//            ret.Dimension = x.Dimension + y.Dimension;
//            ret.Rank = x.Rank + y.Rank - 1;
//            ret.Rank2 = (x.Rank - 1) * (y.Rank - 1) + x.Rank2 + y.Rank2 - 1;
//            if (String.Compare( x.Herkunft,y.Herkunft)>0)
//                ret.Herkunft=x.Herkunft+", "+ y.Herkunft;
//            else
//                ret.Herkunft=y.Herkunft+", "+x.Herkunft;
//            return ret;
//        }
//    }
//}
