using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SqlClient;
using System.Data;
using System.Numerics;
using System.Diagnostics;

namespace EssFinderProject
{
    class CoposPlusCheck
    {
        //private const string connectionString = Settings.DatabaseConnectionString;
        private const int _fkrun = 6607;//6597;
        private const int anzahless = 7;


        public CoposPlusCheck()
        {

            using (SqlConnection conn = new SqlConnection(MainClass.DatabaseConnectionString))
            {
                //conn.ConnectionTimeout = 0;
                conn.Open();
                var cmd = conn.CreateCommand();
                cmd.CommandTimeout = 0;
                //cmd.CommandText = "select matrixid, matrix, dimension, (select top 1 payoff from vektoren v where v.matrixid=m.matrixid and v.fk_run=m.fk_run and isess=1 order by payoffDecimal desc) 'payoff' from matrizen m where anzahless>=" + anzahless + " and m.fk_run =" + _fkrun + " order by matrixid";
                cmd.CommandText = "exec CopositivePlusCheck";
                var rMatizen = cmd.ExecuteReader();
                cmd.Dispose();
                int counter=0;
                while (rMatizen.Read())
                {
                    string line = rMatizen["matrix"].ToString();
                    int dimension = Convert.ToInt32( rMatizen["dimension"]);

                    //MatrixBigRational matbr = MatrixBigRational.MatrixFromCircularDatabaseString(rMatizen["matrix"].ToString());
                    MatrixBigRational matbr = MatrixBigRational.MatrixFromMathematicaString(rMatizen["matrix"].ToString());
                    BigRational gamma = new BigRational(rMatizen["payoff"].ToString());
                    MatrixBigRational x = gamma * MatrixBigRational.OneMatrix(matbr.Rows, matbr.Columns) - matbr;

                    if ((x.SignOfDeterminantMathematica() != 0) && (x.IsPositiveSemidefiniteMathematica() == false))
                    {
                        LinkedList<int> subsetList = BitVector.GetAllSubsetsWithExactelyNElements(dimension, 3);
                        bool detFalsch = false;
                        foreach (int bitmask in subsetList)
                        {
                            MatrixBigRational y = x.GetPrincipalSubmatrixByBitmask(bitmask, 3);
                            if (y.SignOfDeterminantMathematica() == -1)
                            {
                                detFalsch = true;
                                break;
                            }
                        }
                        if (detFalsch == false)
                        {
                            Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+rMatizen["matrix"].ToString());
                            Console.WriteLine(rMatizen["matrixid"]);
                            //Console.WriteLine("Regulär: " + (x.SignOfDeterminantMathematica() != 0));
                            //Console.WriteLine("PosDef: " + x.IsPosivieDefiniteMathematica());
                            //Console.WriteLine("PosSemdef: " + x.IsPositiveSemidefiniteMathematica());
                        }
                    }
                    counter += 1;
                    if (counter % 10000 == 0)
                    {
                        Console.WriteLine(counter);
                        Console.WriteLine(line);
                        
                    }

                }
            }
        }


    }
}
