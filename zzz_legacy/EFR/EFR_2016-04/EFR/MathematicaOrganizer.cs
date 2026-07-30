using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Wolfram.NETLink;

namespace EFR
{
    public class MathematicaOrganizer
    {
        IKernelLink ml;
        private  Object lockObjectAll = new Object();
        public MathematicaOrganizer()
        {
            // This launches the Mathematica kernel:
            ml = MathLinkFactory.CreateKernelLink();
            // Discard the initial InputNamePacket the kernel will send when launched.
            ml.WaitAndDiscardAnswer();
        }

        public string EvaluateToInputForm(string input)
        {
            lock (lockObjectAll)
            {
                string result = ml.EvaluateToInputForm(input, 0);
                return result;
            }
        }

        public string EvaluateToOutputForm(string input)
        {
            lock (lockObjectAll)
            {
                string result = ml.EvaluateToOutputForm(input, 0);
                return result;
            }
        }


        public void Close()
        {
            ml.Close();
        }
    }
}



