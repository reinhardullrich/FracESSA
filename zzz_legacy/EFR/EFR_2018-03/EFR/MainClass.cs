using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.IO;
using System.Threading;
using System.Numerics;


namespace EFR
{
    public static class MainClass
    {
        public static MathematicaOrganizer MO;

        public static string temp =  @"C:\Users\ullrich\_py\efr_in.txt";
        public static string tempm = @"C:\Users\ullrich\_py\efr_out_m.txt";
        public static string tempv = @"C:\Users\ullrich\_py\efr_out_v.txt";
        public static int NumberOfWorkerThreads = 6; //sollte anzahl prozessorkerne sein    
        public static int MaxNumberPCQueue = 100000;
        public static int PCEnqueueWaitingTime = 2000000; //in milli.sec.
        public static bool WithLog = false;

        public static void Run()
        {
            MO = new MathematicaOrganizer();
            EssFinderQueue essQueue = new EssFinderQueue();

            var matrices = System.IO.File.ReadLines(MainClass.temp);
            int matrixid = 1;
            foreach (string strmatrix in matrices)
            {
                if (strmatrix == "")
                    continue;
                string[] splitname = strmatrix.Split(new char[] { ';' });
                string[] splitter = splitname[0].Split(new char[] { '+' });
                int nrow = Convert.ToInt32(splitter[0]);
                string[] splitter2 = (splitter[1]).Split(new char[] { ',' });
                MatrixBigRational x = new MatrixBigRational(nrow, splitter2);
                EssFinder essf = new EssFinder(x);
                essf.MatrixId = matrixid;
                if (splitname.Length==2)
                    essf.Name = splitname[1];
                essQueue.EnqueueItem(essf);
                matrixid++;
            }

            essQueue.Shutdown();
        }
    }
}