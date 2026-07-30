using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EFR
{
    public class zzz_MatrixInformation
    {
        public long MatrixId = -1;
        public string MatrixName = "";
        public string Matrix = "";
        public int Dimension = -1;
        public int AnzahlEss = -1;
        public HashSet<VectorInformation> Vectors = new HashSet<VectorInformation>();
    }
}
