using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EFR
{
    public class VectorInformation
    {
        public long MatrixId = -1;
        public string Vektor = "";
        public string Support = "";
        public int SupportSize = -1;
        public string ExtendedSupport = "";
        public int ExtendedSupportSize = -1;
        public bool IsEss = false;
        public byte GrundEss=32;
        public decimal Payoff = 0;
    }
}
