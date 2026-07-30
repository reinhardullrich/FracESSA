using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace tester
{
    class Program
    {
        static void Main(string[] args)
        {
            var sw = new Stopwatch();
            sw.Start();
            EFR.MainClass.Run();
            //Console.ReadLine();
            sw.Stop();
            Console.WriteLine(sw.Elapsed.TotalSeconds);
        }
    }
}
