using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Data.SqlClient;
using System.Data;

namespace EssFinderProject.RankSumTool
{
    class RankSumToolLight
    {

        private const int MaxDimension=1000;
        private const int MaxIterations = 30;

        private struct RankInfo
        {
            public int Dimension ;
            public int Rank ;
            public int Rank2;
            public string Herkunft;
            public int NumberOfBaseMatrices;
            public int d ;
            public int s ;
            public decimal c;


            public RankInfo( int dimension, int rank, int rank2, string herkunft, int numberOfBaseMatrices, int d, int s)
            {
                this.Dimension = dimension;
                this.Rank = rank;
                this.Rank2 = rank2;
                this.Herkunft = herkunft;
                this.NumberOfBaseMatrices = numberOfBaseMatrices;
                this.d = d;
                this.s = s;
                this.c = 0M;

            }
            
        }

        private RankInfo[,] _RankArray = new RankInfo[MaxIterations+1, MaxDimension];

        public RankSumToolLight()
        {
            Init();

            for (int i = 1; i <= MaxIterations; i++)
            {
                for (int j = 0; j < MaxDimension; j++)
                {
                    _RankArray[i, j] = _RankArray[i - 1, j];
                }
                for (int j = 0; j < MaxDimension; j++)
                    for (int k = 0; k < MaxDimension; k++)
                    {
                        var testRank = Calculate(_RankArray[i - 1, j], _RankArray[0, k]);

                        if (testRank.Dimension<=MaxDimension)
                        {
                            var oldRank= _RankArray[i,testRank.Dimension-1];
                            if (testRank.Rank2 > oldRank.Rank2 || (testRank.Rank2 == oldRank.Rank2 && testRank.Rank > oldRank.Rank))
                                _RankArray[i, testRank.Dimension - 1] = testRank;

                        }
                    }

                //PrintIteration(i);
                Console.WriteLine(i);

            }

            //PrintIteration(MaxIterations);
            for (int i = 0; i < MaxDimension; i++)
                calculateC(ref _RankArray[MaxIterations, i]);

            SaveToDatabase(MaxIterations);
            Console.WriteLine("Fertig!");
        }




        private void Init()
        {
            for (int i = 0; i < MaxDimension;i++ )
            {
                _RankArray[0, i] = new  RankInfo( i + 1, i + 1, i + 1, "I_" + (i + 1),1,d_n(i+1),s_n(i+1));  
            }

            //_RankArray[0, 5] = new RankInfo(6, 6, 8, "U_6", 1, d_n(6), s_n(6));
            //_RankArray[0, 6] = new RankInfo(7, 7, 14, "U_7", 1, d_n(7), s_n(7));
            //_RankArray[0, 7] = new RankInfo(8, 8, 18, "U_8", 1, d_n(8), s_n(8));
            //_RankArray[0, 8] = new RankInfo(9, 9, 26, "U_9", 1, d_n(9), s_n(9));
            //_RankArray[0, 10] = new RankInfo(11, 11, 32, "U_11", 1, d_n(11), s_n(11));
            //_RankArray[0, 12] = new RankInfo(13, 13, 50, "U_13", 1, d_n(13), s_n(13));
            //_RankArray[0, 14] = new RankInfo(15, 15, 95, "U_15", 1, d_n(15), s_n(15));
        }


        private int d_n(int dimension)
        {
            return (int)Math.Floor((decimal)dimension * (decimal)dimension / 4);
        }

        private int s_n(int dimension)
        {
            return dimension * (dimension + 1) / 2 - 4;
        }

        private void calculateC(ref RankInfo rankInfo)
        {
            rankInfo.c = ((rankInfo.s - rankInfo.Rank2 - (decimal)Math.Sqrt(2) * (decimal)Math.Pow(rankInfo.Dimension, 1.5)) / rankInfo.Dimension);
        }


        private RankInfo Calculate(RankInfo x, RankInfo y)
        {
            RankInfo ret = new RankInfo();
            ret.Dimension = x.Dimension + y.Dimension;
            ret.Rank = x.Rank + y.Rank - 1;
            ret.Rank2 = (x.Rank - 1) * (y.Rank - 1) + x.Rank2 + y.Rank2 - 1;
            if (String.Compare( x.Herkunft,y.Herkunft)>0)
                ret.Herkunft=x.Herkunft+", "+ y.Herkunft;
            else
                ret.Herkunft=y.Herkunft+", "+x.Herkunft;
            ret.s = s_n(ret.Dimension);
            ret.d = d_n(ret.Dimension);
            ret.NumberOfBaseMatrices = x.NumberOfBaseMatrices + y.NumberOfBaseMatrices;
            return ret;
        }

        private void PrintIteration(int iteration)
        {
            Debug.WriteLine("Durchlauf: " + iteration);
            for (int i=0;i<MaxDimension;i++)
            {
                var x = _RankArray[iteration, i];
                Debug.WriteLine("Iteration "+iteration +": Dimension: " +x.Dimension + " Rank: " + x.Rank + " Rank2: " + x.Rank2 + " BaseMatrices: " + x.Herkunft + " s_n: " + x.s + " d_n: " + x.d+ (x.Rank2<=x.d? "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!":""));
            }
        }

        private void SaveToDatabase(int iteration)
        {
            SqlBulkCopy bc = new SqlBulkCopy(MainClass.DatabaseConnectionString);
            DataTable dt = new DataTable();
            bc.DestinationTableName = "cpraenge";
            using (SqlConnection conn = new SqlConnection(MainClass.DatabaseConnectionString))
            {
                conn.Open();
                using (var cmd = conn.CreateCommand())
                {
                    cmd.CommandText = "SELECT top 0 [run],[iteration],[dimension],[rank],[rank2],[basematrices],[numberbasematrices],[sn],[dn],[c] FROM [dbo].[cpraenge] ";
                    using (var reader = cmd.ExecuteReader())
                        dt.Load(reader);
                }
            }
            bc.BulkCopyTimeout = 0;
            string yy = DateTime.Now.ToString("yyMMddHmmssfff");
            long runId = Convert.ToInt64(yy);
            for (int i = 0; i < MaxDimension; i++)
            {
                DataRow dr = dt.NewRow();
                var x = _RankArray[iteration, i];
                dr[0] = runId;
                dr[1] = iteration;
                dr[2] = x.Dimension;
                dr[3] = x.Rank;
                dr[4] = x.Rank2;
                dr[5] = x.Herkunft;
                dr[6] = x.NumberOfBaseMatrices;
                dr[7] = x.s;
                dr[8] = x.d;
                dr[9] = x.c;

                dt.Rows.Add(dr);
            }

            bc.WriteToServer(dt);
            dt.Clear();

        }

    }
}
