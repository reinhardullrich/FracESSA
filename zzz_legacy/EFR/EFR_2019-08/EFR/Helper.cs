using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

namespace EFR
{
    public enum ReasonEss : int {

		indeterminate=-1,
        true_pure_ess=1,
        true_posdef_double = 2,
		true_posdef_rational=3,
        true_copositive = 4,
        false_not_posdef_and_K_0_1 = 5,
        false_not_partial_copositive = 6,
        false_not_copositive = 7
        //true_copos_wrt_Kequals0_posdef=4,
        //true_copos_wrt_Kequals2andJequals2=5,

        //false_copos=6,
        //false_copos_wrt_diag_leq_0=7,
        //false_copos_wrt_Kequals2andJequals2 = 8,
        //true_copos_BR=9,
        //false_copos_BR=10,

    }
			
    public static class Helper
    {

        public static string Left(this string str, int count)
        {
            if (string.IsNullOrEmpty(str) || count < 1)
                return string.Empty;
            else
                return str.Substring(0, Math.Min(count, str.Length));
        }
			

        public static void ShiftRight(BigRational[] a)
        {
            BigRational temp = a[a.Length - 1];
            for (int i=a.Length-1;i>0;i--)
                a[i]=a[i-1];
            a[0]=temp;
        }



        public static string ToStringExtended<T>(this IList<T> list)
        {
            return string.Join(",", list.ToArray());
        }

        public static List<int> ListofElement(int i)
        {
            List<int> list = new List<int>();
            list.Add(i);
            return list;
        }
    }
}

//        public static bool IsOdd(int value)
//        {
//            return value % 2 != 0;
//        }



//        public static IEnumerable<T[]> Permutations<T>(T[] values, int fromInd = 0)
//        {
//            if (fromInd + 1 == values.Length)
//                yield return values;
//            else
//            {
//                foreach (var v in Permutations(values, fromInd + 1))
//                    yield return v;
//
//                for (var i = fromInd + 1; i < values.Length; i++)
//                {
//                    SwapValues(values, fromInd, i);
//                    foreach (var v in Permutations(values, fromInd + 1))
//                        yield return v;
//                    SwapValues(values, fromInd, i);
//                }
//            }
//        }
//
//        private static void SwapValues<T>(T[] values, int pos1, int pos2)
//        {
//            if (pos1 != pos2)
//            {
//                T tmp = values[pos1];
//                values[pos1] = values[pos2];
//                values[pos2] = tmp;
//            }
//        }

//        public static BigRational[] ToBRLine(string[] line)
//        {
//            var ret = new BigRational[line.Length];
//            for (int i = 0; i < line.Length; i++)
//                ret[i] = new BigRational(line[i]);
//            return ret;
//        }

//        public static string ToEFName(this int[] input)
//        {
//            StringBuilder sB = new StringBuilder();
//            for (int i = 0; i < input.Length; i++)
//                sB.Append(input[i]).Append(",");
//            sB.Length--;
//            return sB.ToString();
//        }
//
//        public static string ToEFName(this decimal[] input)
//        {
//            StringBuilder sB = new StringBuilder();
//            for (int i = 0; i < input.Length; i++)
//                sB.Append(input[i]).Append(";");
//            sB.Length--;
//            return sB.ToString();
//        }
//
//        public static int[] FromEFName(this string input)
//        {
//            return input.Split(',').Select(s => Convert.ToInt32(s)).ToArray();
//            //return input.Split(new string[] {","}, StringSplitOptions.RemoveEmptyEntries);
//            //StringBuilder sB = new StringBuilder();
//            //for (int i = 0; i < input.Length; i++)
//            //    sB.Append(input[i]).Append(",");
//            //sB.Length--;
//            //return sB.ToString();
//        }
//


//public static decimal SinM(decimal angle)
//{
//    BigRational x = new BigRational(angle);
//    string input = "N[Sin[Pi*" + x.Numerator +"/" + x.Denominator +"],40]";
//    //Debug.Print(input);
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    return Convert.ToDecimal(ret.Replace(".",","));
//}

//public static decimal CosM(decimal angle)
//{
//    BigRational x = new BigRational(angle);
//    string input = "N[Cos[Pi*" + x.Numerator + "/" + x.Denominator + "],40]";
//    //Debug.Print(input);
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    return Convert.ToDecimal(ret.Replace(".", ","));
//}

//public static decimal SinME(int nominator, int denominator, decimal anglemin, decimal anglemax)
//{
//    string input = ("N[Sin[Pi*(" + anglemin.ToString() +"+" + nominator.ToString() + "*("+ anglemax.ToString() +"-" + anglemin.ToString() +")/" + denominator.ToString()).Replace(",",".") + ")],70]";
//    //Debug.Print(input);
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    //Debug.Print(ret);
//    return Convert.ToDecimal(ret.Replace(".", ","));
//}

//public static decimal CosME(int nominator, int denominator, decimal anglemin, decimal anglemax)
//{
//    string input = ("N[Cos[Pi*(" + anglemin.ToString() + "+" + nominator.ToString() + "*(" + anglemax.ToString() + "-" + anglemin.ToString() + ")/" + denominator.ToString()).Replace(",",".") + ")],70]";
//    //Debug.Print(input);
//    string ret = MainClass.MO.EvaluateToOutputForm(input);
//    //Debug.Print(ret);
//    return Convert.ToDecimal(ret.Replace(".", ","));
//}

//public static bool isSubset(List<int> subset, List<int> superset)
//{
//    //return !subset.Except(superset).Any();
//    //var x = subset.Except(superset);
//    //return !x.Any();
//    //return (x.Count == 0);
//    bool isSubset = true;
//    for (int shiftDim = 0; shiftDim < subset.Count; shiftDim++)
//        if (superset.IndexOf(subset[shiftDim]) == -1)
//        {
//            isSubset = false;
//            break;
//        }
//    return isSubset;
//}


//public static LinkedList<List<T>> PowerSet<T>(List<T> startingSet, int minSubsetSize)
//{
//    LinkedList<List<T>> subsetList = new LinkedList<List<T>>();

//    //The set bits of each intermediate value represent unique 
//    //combinations from the startingSet.
//    //We can start checking for combinations at (1<<minSubsetSize)-1 since
//    //_values less than that will not yield large enough subsets.
//    int iLimit = 1 << startingSet.Count;
//    for (int shiftDim = (1 << minSubsetSize) - 1; shiftDim < iLimit; shiftDim++)
//    {
//        //Get the number of 1's in this 'shiftDim'
//        int setBitCount = NumberOfSetBits(shiftDim);

//        //Only include this subset if it will have at least minSubsetSize members.
//        if (setBitCount >= minSubsetSize)
//        {
//            List<T> subset = new List<T>(setBitCount);

//            for (int j = 0; j < startingSet.Count; j++)
//            {
//                //If the j'th bit in shiftDim is set, 
//                //then add the j'th element of the startingSet to this subset.
//                if ((shiftDim & (1 << j)) != 0)
//                {
//                    subset.Add(startingSet[j]);
//                }
//            }
//            subsetList.AddLast(subset);
//        }
//    }
//    return subsetList;
//}

//private static int NumberOfSetBits(int shiftDim)
//{
//    shiftDim = shiftDim - ((shiftDim >> 1) & 0x55555555);
//    shiftDim = (shiftDim & 0x33333333) + ((shiftDim >> 2) & 0x33333333);
//    return (((shiftDim + (shiftDim >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
//}
//public static string ToStringLL(this LinkedList<int> list)
//{
//    string x = "";
//    foreach (int l in list)
//        x += BitVector.ToListOfInt(l).ToStringExtended()+"\r\n";
//    x=x.Left(x.Length);
//    return x;
//}

//public static string ToStringLL(this List<Vector> list)
//{
//    string x = "";
//    foreach (Vector l in list)
//        x += l.ToString() + "\r\n";
//    x = x.Left(x.Length);
//    return x;
//}