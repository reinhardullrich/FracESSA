using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace EFR
{
    class BitVector
    {
        public static List<int>[] GetAllSubsetsWithMoreThanZeroElementsArrayList(int dimension)
        {
            List<int>[] array = new List<int>[dimension];
            for (int i = 0; i < dimension; i++)
                array[i] = new List<int>();
            for (int i = 1; i < (1 << dimension); i++)
                array[SupportSize(i) - 1].Add(i);
            return array;
        }

        public static List<int> ToListOfInt(int bitMask, int dimension =32)
        {
            List<int> supportList = new List<int>();
            for (int i=0;i<dimension;i++)
                if ((bitMask & (1<<i))!=0)
                    supportList.Add(i+1);
            return supportList;
        }

        public static string ToListString(int bitMask, int dimension = 32)
        {
            StringBuilder output = new StringBuilder();
            for (int i = 0; i < dimension; i++)
                if ((bitMask & (1 << i)) != 0)
                    output.Append(i + 1).Append(",");
            output.Length--;
            return output.ToString();
        }

        public static int GetSmallestSetFlag(int bitMask, int dimension=32)
        {
            for (int i = 0; i < dimension; i++)
                if ((bitMask & (1 << i)) != 0)
					return i ;
            return -1;
        }

        public static int SupportSize(int i)
        {
            i = i - ((i >> 1) & 0x55555555);
            i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
            return ((((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
        }			
    }
}

//        public static LinkedList<int> GetAllSubsetsWithMoreThanZeroElements(int dimension)
//        {
//            LinkedList<int> subsetList = new LinkedList<int>();
//            for (int i = 1; i < (1 << dimension); i++)
//                subsetList.AddLast(i);
//            return subsetList;
//        }

//public static LinkedList<int>[] GetAllSubsetsWithMoreThanOneElementAndWithOnlyThisSupportArray(int dimension, int support)
//{
//    LinkedList<int>[] array = new LinkedList<int>[dimension + 1];

//    for (int i = 2; i <= dimension; i++)
//        array[i] = new LinkedList<int>();

//    for (int i = 3; i < (1 << dimension); i++)
//        if (((i & (i - 1)) != 0) && ((i & support) == i))
//        {
//            array[SupportSize(i)].AddLast(i);
//        }
//    return array;
//}

//public static LinkedList<int>[] GetAllSubsetsWithMoreThanZeroElementsArray(int dimension)
//{
//    LinkedList<int>[] array = new LinkedList<int>[dimension];
//    for (int i = 0; i < dimension; i++)
//        array[i] = new LinkedList<int>();
//    for (int i = 1; i < (1 << dimension); i++)
//        array[SupportSize(i)-1].AddLast(i);
//    return array;
//}

//public static int ShiftedLeft(int bitMask, int dimension)
//{
//    bool lowestBit = (bitMask & (1 << (dimension - 1))) != 0;
//    bitMask=bitMask << 1;
//    if (lowestBit)
//        bitMask = bitMask | 1;
//    bitMask = bitMask & (~(1 << dimension));
//    return bitMask;
//}


//public static bool IsSubset(int setToCheck, int theSet)
//{
//    return ((setToCheck & theSet) == setToCheck);
//}

//public static int AddElement(int bitMask, int element)
//{
//    return (bitMask | (1 << (element - 1)));
//}



//public static int FromListOfInt(List<int> support)
//{
//    int bitmask = 0;
//    foreach (int supportPosition in support)
//    {
//        bitmask = bitmask | (1 << (supportPosition - 1));
//    }
//    //Debug.WriteLine(Convert.ToString(bitmask,2));
//    return bitmask;
//}

//public static LinkedList<int> GetAllSubsetsWithMoreThanOneElementAndWithOnlyThisSupport(int dimension,int support)
//{
//    LinkedList<int> subsetList = new LinkedList<int>();
//    for (int i = 3; i < (1 << dimension); i++)
//        if (((i & (i - 1)) != 0) && ((i&support)==i))
//        {
//            subsetList.AddLast(i);
//            //Debug.WriteLine(Convert.ToString(shiftDim, 2));
//        }
//    return subsetList;
//}

//public static HashSet<int> GetAllSubsetsWithMoreThanOneElementAndWithOnlyThisSupportHashSet(int dimension, int support)
//{
//    HashSet<int> subsetList = new HashSet<int>();
//    for (int i = 3; i < (1 << dimension); i++)
//        if (((i & (i - 1)) != 0) && ((i & support) == i))
//        {
//            subsetList.Add(i);
//            //Debug.WriteLine(Convert.ToString(shiftDim, 2));
//        }
//    return subsetList;
//}

//public static HashSet<int> GetAllSubsetsWithMoreThanOneElementAndWithOnlyThisSupportLookup(int dimension, int support)
//{
//    Lookup<int,int> subsetList = new Lookup<int, int>();
//    for (int shiftDim = 3; shiftDim < (1 << dimension); shiftDim++)
//        if (((shiftDim & (shiftDim - 1)) != 0) && ((shiftDim & support) == shiftDim))
//        {
//            subsetList.Add(shiftDim);
//            //Debug.WriteLine(Convert.ToString(shiftDim, 2));
//        }
//    return subsetList;
//}

//public static LinkedList<int> GetAllSubsetsWithMoreThanOneElement(int dimension)
//{
//    LinkedList<int> subsetList = new LinkedList<int>();
//    for (int i = 3; i < (1 << dimension); i++)
//        if ((i & (i - 1)) != 0)
//        {
//            subsetList.AddLast(i);
//            //Debug.WriteLine(Convert.ToString(shiftDim, 2));
//        }
//    return subsetList;
//}

//public static LinkedList<int> GetAllSubsetsWithExactelyNElements(int dimension, int n)
//{
//    LinkedList<int> subsetList = new LinkedList<int>();
//    for (int i = 1; i < (1 << dimension); i++)
//        if (NumberOfSetBits(i)==n)
//        {
//            subsetList.AddLast(i);
//            //Debug.WriteLine(Convert.ToString(shiftDim, 2));
//        }
//    return subsetList;
//}