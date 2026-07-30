using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EFR
{
    public class VectorInformation
    {
        public long VectorID;
        public string Vector;
        public string Support;
        public int SupportSize;
        public string ExtendedSupport;
        public int ExtendedSupportSize;
        public int IsEss;
        public int ReasonEss;
		public BigRational Payoff;
        public decimal PayoffDecimal;

		public string GetVectorInformation()
		{
			StringBuilder str = new StringBuilder("");
			str.Append("'" + this.VectorID + "',"); 
			str.Append("'" + this.Vector + "',"); 
			str.Append("'" + this.Support + "',");
			str.Append("'" + this.SupportSize + "',");
			str.Append("'" + this.ExtendedSupport + "',");
			str.Append("'" + this.ExtendedSupportSize + "',");
			str.Append("'" + this.IsEss + "',");
			str.Append("'" + this.ReasonEss + "',");
			str.Append("'" + this.Payoff.ToString() + "',");
			str.Append("'" + this.PayoffDecimal.ToString().Replace(",", ".") + "'");
			return(str.ToString ());

		}


		public static string GetHeader()
		{
			StringBuilder str = new StringBuilder ("");
			str.Append ("'VectorID',");
			str.Append ("'Vector',");
			str.Append ("'Support',");
			str.Append ("'SupportSize',");
			str.Append ("'ExtendedSupport',");
			str.Append ("'ExtendedSupportSize',");
			str.Append ("'IsEss',");
			str.Append ("'Reason',");
			str.Append ("'Payoff',");
			str.Append ("'PayoffDecimal'");
			return(str.ToString ());
		}
    }
}
