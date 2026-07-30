#include <iostream>
#include <vector>
#include <cstdint>
#include <cassert>
#include <string>


#include "include/EssFinder.h"
#include "include/MatrixTemplate.h"
#include "argtable3.h"



struct arg_lit *vectors, *help, *version, *logger, *exact, *fullsupport;
struct arg_str *matrixinput;
struct arg_end *end;

bool parseMatrix(std::string input, MatrixTemplate<Rational>& matrix);

int main(int argc, char *argv[])
{
//    argc = 3;
//    //argv[1]="3#1,2,3,7,1,7/2,9,11/3,0";
//    //argv[1]="7#2/3,5,9";
//    //argv[1]="9#2,2,5,9";
//    //argv[1]="21#15,15,7,15,15,7,7,15,7,15";
//    //argv[1]="23#27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939";
//    //argv[1]="24#15,15,7,15,15,7,15,7,7,15,15,7";
//    //argv[1]="15#9,9,5,9,3,5,9";
//    //argv[1]="22#2,2,4,4,5,5,4,4,2,2,-1";
//    //argv[1]="4#-1,1,0,0,1,-1,0,0,0,0,0,-1,0,0,-1,0";
//    //argv[1]="3#0,0,0,0,0,0,1,1,0";
//    //argv[1]="2#1,2,3,4";
//    argv[1]="5#1,0,2,2,2,0,1,2,2,2,2,2,1,0,0,2,2,0,1,0,2,2,0,0,0";
//    //argv[1]="19#1,2,2,2,2,2,1,1,2";
//    //argv[1]="30#9,9,5,9,3,5,9,9,3,9,5,9,3,5,9";
//    //argv[1]="12#287,50,200,150,113,200";
//    //argv[1]="12#292,-1,293,-1,293,-1";
//    //argv[1]="12#16029525,4007380,9088537,12022143,2147548,16029525";
//    //argv[1]="19#13,31,31,27,31,27,13,13,27";
//    //argv[1]="5#1,3";
//    argv[2]= "-v";


    /* the global arg_xxx structs are initialised within the argtable */
    void *argtable[] = {
        help            = arg_litn(NULL, "help", 0, 1, "display this help and exit"),
        version         = arg_litn(NULL, "version", 0, 1, "display version info and exit"),
        vectors         = arg_litn("v", "vectors", 0, 1, "include the resulting vectors in the output"),
        logger          = arg_litn("l", "log", 0, 1, "output a detailed log file named 'log.txt' in the directory of the program, for diagnostic of learning purposes only"),
        exact           = arg_litn("e", "exact", 0, 1, "only uses rational numbers, for matrices with extreme differences in the input, is much much slower!"),
        fullsupport     = arg_litn("f", "fullsupport", 0, 1, "searches the full support directly after searching support size one. Enable if you expect the matrix to have exactly one ess in the interior of the simplex!"),
        matrixinput     = arg_strn(NULL, NULL, "<matrix>", 1, 1, "the matrix to compute"),
        end     = arg_end(50),
    };

    int exitcode = 0;
    char progname[] = "REF.exe";

    int nerrors;
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        printf("Demonstrate command-line parsing in argtable3.\n\n");
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        exitcode = 0;
        goto exit;
    }
    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout, end, progname);
        printf("Try '%s --help' for more information.\n", progname);
        exitcode = 1;
        goto exit;
    }

    if (nerrors == 0)
    {
        config::withLog=logger->count;
        config::withVectors=vectors->count;
        config::exact=exact->count;
        config::fullSupport=fullsupport->count;

        MatrixTemplate<Rational> A;
        exitcode=parseMatrix(matrixinput->sval[0], A);

        if (exitcode==0) {
            EssFinder x = EssFinder(A);
            x.computeEss();

//            std::cout<<smallestRepresentation(11,6);
//            std::cout<<smallestRepresentation(43,6);
//            Rational myints[] = {-14423932725568825833984914012016456550789939200000,28847865451137651667969828024032913101579878400000,-14423932725568825833984914012016456550789939200000,-14423932725568825833984914012016456550789939200000,28847865451137651667969828024032913101579878400000,28847865451137651667969828024032913101579878400000};
//            std::vector<Rational> y (myints, myints + sizeof(myints) / sizeof(Rational) );
//            MatrixTemplate<Rational> A  = MatrixTemplate<Rational>(y,1);
//            A.printMatrix();
//            std::cout<<A.determinant()<<std::endl;//.printMatrix();
//            std::cout<<A.determinantLaplace()<<std::endl;
//            std::cout<<A.isCopositive();//.adjugate().printMatrix();
//            std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
//            std::cout<<A.isCopositiveLaplace();//.adjugateLaplace().printMatrix();
//            A.multiplyWith(A.adjugate()).printMatrix();
//            std::cout<<A.determinant()<<std::endl;//.printMatrix();
//            std::cout<<A.determinantLaplace();
            //std::cout<<A.toDouble().determinant();
            //MatrixTemplate<double> xx = MatrixTemplate.IdentityMatrix(5);
        }
        goto exit;
    }

exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return exitcode;
}

bool parseMatrix(std::string input, MatrixTemplate<Rational>& matrix)
{
    if (input.find("#")==std::string::npos) {
        std::cout<<"ERROR: string for the matrix does not include '#' as a seperator between dimension and matrix!"<<std::endl;
        return 1;
    }
    std::vector<std::string> dimensionMatrix=splitString(input,'#');
    size_t dimension = std::stoull(dimensionMatrix[0]);

    std::vector<std::string> matrixStringArray = splitString(dimensionMatrix[1],',');

    std::vector<Rational> matrixArray;

    for (auto x: matrixStringArray) {
        matrixArray.push_back(Rational(x));
    }

    if (!(matrixArray.size()==dimension/2 || matrixArray.size()==dimension*dimension)) {
        std::cout<<"ERROR: The number of matrix-elements must either be floor(dimension/2) (for a cyclically symmetric matrix) or dimension^2 (for any other matrix)! Neither of that is the case!"<<std::endl;
        return 1;
    }

    if (matrixArray.size()==dimension/2) { //cyclically symmetric
        config::isCS = true;
        matrix=MatrixTemplate<Rational> (dimension, matrixArray);
    } else
        matrix=MatrixTemplate<Rational> (matrixArray);
    return 0;
}
