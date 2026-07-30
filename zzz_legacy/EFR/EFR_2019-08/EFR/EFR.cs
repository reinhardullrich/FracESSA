using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
//using NDesk.Options;
using System.IO;
using CommandLine;
using CommandLine.Text;

namespace EFR
{
	class EfrOptions
	{
		[Option('l', null,  HelpText = "Output detailed log-file named in 'Log.txt' in the folder where the .exe is executed.")]
		public bool WithLog { get; set; }

		[Option('v', null, HelpText = "Includes the vectors in a .csv-like manner starting at the second line of the output.")]
		public bool IncludeVectors { get; set; }

		[ValueOption(0)]
		public string MatrixString { get; set; }
	}


    class EFR
    {
		public static bool WithLog = true;

        static void Main(string[] args)
        {
			//args = new string[] { "--help" };
			//args= new string[] {"-lv","4#-1,1,0,5,1,-1,5,0,3,4,0,-1,0,0,-1,0"};
			//args= new string[] {"-v","4#-1,1,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0"};
			//            var sw = new Stopwatch();
			//            sw.Start();

			var options = new EfrOptions();

			if (!CommandLine.Parser.Default.ParseArguments (args, options)) {
				Console.WriteLine ("ERROR: You got an error with your command line options! Please check http://...");
				return;
			}

			EFR.WithLog=options.WithLog;

			if (options.MatrixString == null) {
				Console.WriteLine ("ERROR: No matrix supplied!");
				return;
			}

			MatrixBigRational x=null;
			if (options.MatrixString.IndexOf ("#") > 0) {
				string[] splitter = options.MatrixString.Split(new char[] { '#' });
				string[] matrixarray = (splitter [1]).Split (new char[] { ',' });
				int dimension = Convert.ToInt32 (splitter [0]);
				if (matrixarray.Length==dimension/2)
					x = new MatrixBigRational(dimension,matrixarray);
				else if (matrixarray.Length==dimension*dimension)
					x = new MatrixBigRational(matrixarray);
				else {
					Console.WriteLine("ERROR: Either your dimension our the number of matrix-elements is wrong!");
					return;
				}
			} else {
				Console.WriteLine ("ERROR: The structur of your matrix input does not satisfy: [dimension]#[matrixelements]!");
			}

			EssFinder ef = new EssFinder(x);
			ef.FindESS();
            //Console.BufferHeight = Int16.MaxValue - 1;
			Console.WriteLine (ef.NumberOfEss);
			if (options.IncludeVectors)
				Console.WriteLine (ef.GetVectorCsv ());
			Console.Out.Flush();

			if (EFR.WithLog) {
				string path = System.IO.Path.GetDirectoryName (System.Reflection.Assembly.GetExecutingAssembly ().Location);
				StreamWriter logwriter = new StreamWriter (path + "/Log.txt");
				logwriter.Write (ef.Log);
				logwriter.Flush ();
				logwriter.Close ();
				logwriter.Dispose ();
			}
				
		}
    }
}

//		public static string temp = Environment.CurrentDirectory + "/efr_in.txt"; //@"C:\Users\ullrich\_py\efr_in.txt";
//		public static string tempm = Environment.CurrentDirectory + "/efr_out_m.txt"; //@"C:\Users\ullrich\_py\efr_out_m.txt";
//		public static string tempv = Environment.CurrentDirectory + "/efr_out_v.txt"; //@"C:\Users\ullrich\_py\efr_out_v.txt";
//		public static int NumberOfWorkerThreads = 1; //sollte anzahl prozessorkerne sein    
//		public static int MaxNumberPCQueue = 100000;
//		public static int PCEnqueueWaitingTime = 2000000; //in milli.sec.


//			var matrices = System.IO.File.ReadLines(temp);
//			int matrixid = 1;
//			foreach (string strmatrix in matrices)
//			{
//				if (strmatrix == "")
//					continue;
//				string[] splitname = strmatrix.Split(new char[] { ';' });
//				string[] splitter = splitname[0].Split(new char[] { '+' });
//				int nrow = Convert.ToInt32(splitter[0]);
//				string[] splitter2 = (splitter[1]).Split(new char[] { ',' });
//				MatrixBigRational x = new MatrixBigRational(nrow, splitter2);
//				EssFinder essf = new EssFinder(x);
//				essf.MatrixId = matrixid;
//				if (splitname.Length==2)
//					essf.Name = splitname[1];
//				essQueue.EnqueueItem(essf);
//				matrixid++;

//			EssFinder essf = new EssFinder (x);
//			essQueue.EnqueueItem (essf);
//			essQueue.Shutdown();
//			}
