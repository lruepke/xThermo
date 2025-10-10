#ifndef SWEOSBASH_H
#define SWEOSBASH_H
#include "config.h"
#include "getopt_arguments.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

#include "H2ONaCl.h"
#include "MultiProgressBar.h"
#ifndef USE_OMP
#else
    #include "omp.h"
#endif

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

// #ifdef _WIN32
//     #include "windows.h"
//     #define BLACK			0
//     #define BLUE			1
//     #define GREEN			2
//     #define CYAN			3
//     #define RED				4
//     #define MAGENTA			5
//     #define BROWN			6
//     #define LIGHTGRAY		7
//     #define DARKGRAY		8
//     #define LIGHTBLUE		9
//     #define LIGHTGREEN		10
//     #define LIGHTCYAN		11
//     #define LIGHTRED		12
//     #define LIGHTMAGENTA	13
//     #define YELLOW			14
//     #define WHITE			15
//     static HANDLE   m_hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
//     static WORD     m_currentConsoleAttr;
//     static CONSOLE_SCREEN_BUFFER_INFO csbi;
//     #define COLOR_PURPLE ""
//     #define COLOR_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<""
//     #define COLOR_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (GREEN & 0x0F) );cout<<""
//     #define COLOR_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<""
//     #define COLOR_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BLUE & 0x0F) );cout<<""
//     #define COLOR_DEFAULT "";SetConsoleTextAttribute(m_hConsole, m_currentConsoleAttr );cout<<""
//     #define ERROR_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<"Error: "<<COLOR_DEFAULT
//     #define WARN_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<"Warning: "<<COLOR_DEFAULT
// #else
//     // define color, this seems only work on MacOS and linux, doesn't work on windows
//     #define ERROR_COUT "["<<"\033[31mError: "<<"\033[0m] "
//     #define WARN_COUT "["<<"\033[33mWarning: "<<"\033[0m] "
//     #define COLOR_PURPLE "\033[35m"
//     #define COLOR_RED "\033[31m"
//     #define COLOR_GREEN "\033[32m"
//     #define COLOR_YELLOW "\033[33m"
//     #define COLOR_BLUE "\033[34m"
//     #define COLOR_DEFAULT "\033[0m"
// #endif

// ============= Constants ====================
#ifndef Kelvin
#define Kelvin 273.15
#endif
// ============================================


namespace SWEOSbash
{                                                                                            
    bool bash_run(int argc, char** argv);
    // calculation mode: 0d, 1d, 2d, 3d
    #define CALCULATION_MODE_SINGLEPOINT 0
    #define CALCULATION_MODE_ONEDIMENSION 1
    #define CALCULATION_MODE_TWODIMENSION 2
    #define CALCULATION_MODE_THREEDIMENSION 3
    
    // variable selection
    #define VARIABLE_SELECTION_PTX 0
    #define VARIABLE_SELECTION_PHX 1
    #define VARIABLE_SELECTION_T 2
    #define VARIABLE_SELECTION_P 3
    #define VARIABLE_SELECTION_X 4
    #define VARIABLE_SELECTION_H 5
    #define VARIABLE_SELECTION_PT 6
    #define VARIABLE_SELECTION_PX 7
    #define VARIABLE_SELECTION_TX 8
    #define VARIABLE_SELECTION_PH 9
    #define VARIABLE_SELECTION_HX 10

    #define PROGRAM_THERMO 1
    #define PROGRAM_LUTGEN 2
    #define PROGRAM_LUTINFO 3
    #define PROGRAM_LUT2VTU 4
    #define PROGRAM_LUT2PI 5
    #define PROGRAM_LOOKUP 6

    #define FLUID_H2O 1
    #define FLUID_NaCl 2
    #define FLUID_H2ONaCl 3

    #define BACKEND_H2O_IAPS84 1
    #define BACKEND_H2O_IAPWS95 2
    #define BACKEND_H2O_IAPWS95_COOLPROP 3

    static void StartText()
    {
        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<COLOR_YELLOW;       //print text in yellow color
        cout << "***************************************************\n";
        cout << "*                 program xThermal                   *\n";
        cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
        cout << "*  Version: "<<xThermal_VERSION<<"                     *\n";
        cout << "*                                                 *\n";
        cout << "*  Equation of state of H2O, NaCl, H2O-NaCl     *\n";
        cout << "*  - Independent variables: TPX, HPX              *\n";
        cout << "*  - Properties: density, enthalpy, viscosity     *\n";
        cout << "*  - saturation, salinity, phase diagram          *\n";
        cout << "*  unit:                                          *\n";
        cout << "*      temperature-deg.C,        pressure-Pa     *\n";
        cout << "*      salinity-kg/kg,       density-kg/m3       *\n";
        cout << "*      enthalpy-J/kg,        viscosity-Pa s      *\n";
        cout << "*                                                 *\n";
        cout << "* (c) Zhikui Guo, GEOMAR, "<<xThermal_DATE<<", Kiel        *\n";
        cout << "*                                                 *\n";
        cout << "***************************************************\n";
        cout << "\n";
        cout<<COLOR_DEFAULT;

    }

    static void StartText_artASCII()
    {
        // cout<<COLOR_GREEN<<"███████╗ █████╗ ██╗  ████████╗██╗    ██╗ █████╗ ████████╗███████╗██████╗     ███████╗ ██████╗ ███████╗\n"
        // <<"██╔════╝██╔══██╗██║  ╚══██╔══╝██║    ██║██╔══██╗╚══██╔══╝██╔════╝██╔══██╗    ██╔════╝██╔═══██╗██╔════╝\n"
        // <<"███████╗███████║██║     ██║   ██║ █╗ ██║███████║   ██║   █████╗  ██████╔╝    █████╗  ██║   ██║███████╗\n"
        // <<"╚════██║██╔══██║██║     ██║   ██║███╗██║██╔══██║   ██║   ██╔══╝  ██╔══██╗    ██╔══╝  ██║   ██║╚════██║\n"
        // <<"███████║██║  ██║███████╗██║   ╚███╔███╔╝██║  ██║   ██║   ███████╗██║  ██║    ███████╗╚██████╔╝███████║\n"
        // <<"╚══════╝╚═╝  ╚═╝╚══════╝╚═╝    ╚══╝╚══╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝    ╚══════╝ ╚═════╝ ╚══════╝\n"
        // <<COLOR_DEFAULT<<std::endl;;

        // cout<<COLOR_GREEN<<"                              $$$$$$$$\\  $$$$$$\\   $$$$$$\\  \n"
        // <<"                              $$  _____|$$  __$$\\ $$  __$$\\ \n"
        // <<" $$$$$$$\\ $$\\  $$\\  $$\\       $$ |      $$ /  $$ |$$ /  \\__|\n"
        // <<"$$  _____|$$ | $$ | $$ |      $$$$$\\    $$ |  $$ |\\$$$$$$\\  \n"
        // <<"\\$$$$$$\\  $$ | $$ | $$ |      $$  __|   $$ |  $$ | \\____$$\\ \n"
        // <<" \\____$$\\ $$ | $$ | $$ |      $$ |      $$ |  $$ |$$\\   $$ |\n"
        // <<"$$$$$$$  |\\$$$$$\\$$$$  |      $$$$$$$$\\  $$$$$$  |\\$$$$$$  |\n"
        // <<"\\_______/  \\_____\\____/       \\________| \\______/  \\______/ \n"
        // <<COLOR_DEFAULT<<std::endl;;
    }

    class cSWEOSarg
    { 
    private:
        xThermal::cxThermal* m_pEOS;
        string m_programName,m_fluidName;
        string m_fileName_LUT;//need by lookup program
        map<string, int> m_map_programName2Index, m_map_fluidName2Ind, m_map_backendName2Index;
        bool m_haveD, m_haveV, m_haveP, m_haveT, m_havet, m_haveX, m_haveH, m_haveR, m_haveG, m_haveO, m_haveB, m_haveF, m_haveL, m_havep;
        int m_valueD, m_threadNumOMP;
        string m_valueV, m_valueG, m_valueO, m_valueB;
        double m_valueT, m_valueP, m_valueX, m_valueH;
        bool m_normalize_vtk;
        int m_min_level, m_max_level; //for lutgen
        LOOKUPTABLE_FOREST::EOS_ENERGY m_TorH;
        int m_valueL; //minimum level interface
        // min/delta/max, order coresponding to -V parameter, 
        //e.g. -VPT, m_valueR1 for pressure, m_valueR2 for temperature
        // double m_valueR1[3], m_valueR2[3], m_valueR3[3];
        double m_valueR[3][3];
        vector<string> m_valueR_str;
    public:
        cSWEOSarg(/* args */);
        ~cSWEOSarg();
    public:
        int m_CalculationMode;
        int m_VariableSelection;
    public:
        bool Parse(int argc, char** argv); //Parse arguments
        bool CheckRange_T(double T_C);
        bool CheckRanges_T(double Trange[2]);
        bool CheckRange_P(double P0);
        static bool CheckRange_P(double P0, double Pmin, double Pmax, const std::string& info="");
        static bool CheckRange_TorH(LOOKUPTABLE_FOREST::EOS_ENERGY TorH,double TorH0, double TorHmin, double TorHmax, const std::string& info);
        bool CheckEOS_Energy();
        bool CheckRanges_P(double Prange[2]);
        static bool CheckRange_X(double X0, double XMIN=xThermal::H2ONaCl::X_MIN, double XMAX=xThermal::H2ONaCl::X_MAX, const std::string& info="");
        bool CheckRanges_X(double Xrange[2], double XMIN=0, double XMAX=1);
        bool CheckRange_H(double H0, double P0, double X0);//PHX 0D calculation
        bool CheckRanges_H(double Hrange[2], double P0, double X0);//PHX 0D calculation
        bool CheckRanges_H_PX(double HMIN0, double HMAX0, double PXrange[4]);//PHX 3D calculation
        bool CheckRanges_H_P(double HMIN0, double HMAX0, double Prange[2], double X);//PH and fixed X: 2D calculation
        bool CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P);//HX and fixed P: 2D calculation

        bool checkFluid(string fluidName);
        bool checkProgram(string programName);
        bool Validate(int argc, char** argv); // validate arguments and print corresponding error information
        bool Validate_0D();
        bool Validate_1D();
        bool Validate_2D();
        bool Validate_3D();
        bool Validate_lutinfo(int argc, char** argv);
        bool Validate_lut2vtu(int argc, char** argv);
        bool Validate_lut2pi(int argc, char** argv);
        xThermal::ThermodynamicProperties calculateSinglePoint_PTX(double P, double T_K, double X, bool isCout=true);
        xThermal::ThermodynamicProperties calculateSinglePoint_PHX(double P, double H, double X, bool isCout=true);
        xThermal::ThermodynamicPropertiesVector calculateMultiPoints_PTX_PHX(const string& valueV, const string& filePTX, string outFile, const string& isT_H);
        xThermal::PhaseRegion lookup_TorHPX(const int& dim, const double& TorH, const double& P, const double& X, double* props, bool isCout= false);
        std::vector<xThermal::propInfo> lookup_TorHPX(const std::vector<double>& TorH, const std::vector<double>& P, const std::vector<double>& X, std::vector<std::vector<double>>& props, std::string filename_out_csv="", bool isMeshGrid=false);
        static void helpINFO_thermo()
        {
            STATUS("Help information of "<<COLOR_PURPLE<<"thermo"<<COLOR_GREEN<<" program");
            STATUS("Calculate thermodynamic properties and phase regions for specific input fluid and given TPX or HPX input.");
            STATUS("[Example]: xThermal thermo -D2 -VXP -R0.001/0.001/1/1E5/5E5/2200E5 -T600");
            int wordWidth=20;
            cout<<"options:"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -h "<<COLOR_BLUE<<"List descriptions of usage and available arguments"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -v "<<COLOR_BLUE<<"Print xThermal version number"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -D "<<COLOR_BLUE<<"Dimension: 0, 1, 2, 3. e.g.: -D2"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -V "<<COLOR_BLUE<<"Select independent variables according to -D arguments."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Combination of: T, P, X, H. e.g.: -VXT for 2D and -VTPX for 3D"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -T "<<COLOR_BLUE<<"Set fixed temperature value if T is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -P "<<COLOR_BLUE<<"Set fixed pressure value if P is not in -V option. -P316"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -X "<<COLOR_BLUE<<"Set fixed salinity value if X is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -H "<<COLOR_BLUE<<"Set fixed enthalpy value if H is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -R "<<COLOR_BLUE<<"Set range and interval of variables in -V option, must in the save order with -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"e.g.: -R0.001/0.001/1/0.2/1/1000 for -VXT"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -G "<<COLOR_BLUE<<"Set input filename of TPX or HPX text file for multi-points calculation"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"only used when -D0 and no -P, -X, -T or -H arguments."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"The text file with three columns, PTX or PHX are decided by -V options."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -O "<<COLOR_BLUE<<"Set out put file name, file format is determined by file extension name."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Supported file format is vtk, csv, txt."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -n "<<COLOR_BLUE<<"If normalize the result in vtk file. Used with -D2 or -D3"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -t "<<COLOR_BLUE<<"Set number of thread for parallel computing."<<COLOR_DEFAULT<<std::endl;
            cout<<"Units:"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Temperature "<<COLOR_BLUE<<"Degree Celsius: 1 deg.C = 273.15 K (Kelvin)"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Pressure "<<COLOR_BLUE<<"Pa"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Salinity "<<COLOR_BLUE<<"mass fraction in range [0,1]: seawater is 0.032 = 3.2 wt. % NaCl"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Enthalpy "<<COLOR_BLUE<<"SI: J/kg"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Density "<<COLOR_BLUE<<"SI: kg/m3"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Viscosity "<<COLOR_BLUE<<"SI: Pa s"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  Saturation "<<COLOR_BLUE<<"mass fractionin range of [0, 1]"<<COLOR_DEFAULT<<std::endl;
        }
        void helpINFO_lookup()
        {
            STATUS("Help information of "<<COLOR_PURPLE<<"lookup"<<COLOR_GREEN<<" program");
            STATUS("Lookup thermodynamic properties of input TPX or HPX point from AMR-LUT.");
            STATUS("[Example]: xThermal lookup lut_constT_XP_10.bin -VXP -R0.001/0.001/1/1E5/5E5/2200E5");
            int wordWidth=20;
            cout<<"options:"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -h "<<COLOR_BLUE<<"List descriptions of usage and available arguments"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -v "<<COLOR_BLUE<<"Print xThermal version number"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -D "<<COLOR_BLUE<<"Dimension: 0, 1, 2, 3. e.g.: -D2"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -V "<<COLOR_BLUE<<"Set order of independent variables for multiple points input from -G or input range from -R."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Combination of: T, P, X, H. e.g., -VTPX, three columns of the input file by -G must in order of T-P-X. or -RTmin/dT/Tmax/Pmin/dP/Pmax/Xmin/dX/Xmax"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -T "<<COLOR_BLUE<<"Set fixed temperature value if T is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -P "<<COLOR_BLUE<<"Set fixed pressure value if P is not in -V option. -P316"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -X "<<COLOR_BLUE<<"Set fixed salinity value if X is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -H "<<COLOR_BLUE<<"Set fixed enthalpy value if H is not in -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -R "<<COLOR_BLUE<<"Set range and interval of variables in -V option, must in the save order with -V option."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"e.g.: -R0.001/0.001/1/0.2/1/1000 for -VXT"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -G "<<COLOR_BLUE<<"Set input filename of TPX or HPX text file for multi-points calculation"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"only used when -D0 and no -P, -X, -T or -H arguments."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"The text file with three columns, PTX or PHX are decided by -V options."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -O "<<COLOR_BLUE<<"Set out put file name, file format is determined by file extension name."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Supported file format is vtk, csv, txt. Note that the vtk format is only valid for mesh grid calculation (-R option)."<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -n "<<COLOR_BLUE<<"If normalize the result in vtk file. Used with -R option"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -t "<<COLOR_BLUE<<"Set number of thread for parallel computing."<<COLOR_DEFAULT<<std::endl;
        }
        void helpINFO_Program()
        {
            STATUS("Help information of xThermal main program");
            string version=xThermal_VERSION;
            string author="Zhikui Guo";
            string locus="GEOMAR, Germany";
            string email="zguo@geomar.de";
            int wordWidth=20;
            // time_t now=time(0);
            // char* now_str=ctime(&now);
            string now_str=xThermal_DATE;

            //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
            cout<<"========================== xThermal ==========================="<<std::endl;;
            cout<<"xThermal, a multi-platform program for x-fluid Equation of State and thermodynamic properties calculation. x=[H2O | H2O-NaCl]"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author "<<COLOR_GREEN<<author<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus "<<COLOR_GREEN<<locus<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date "<<COLOR_GREEN<<now_str<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version "<<COLOR_GREEN<<version<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Email "<<COLOR_GREEN<<email<<COLOR_DEFAULT<<std::endl;;
            cout<<"============================================================"<<std::endl;;
            cout<<COLOR_BLUE<<"Usage:   xThermal ["<<COLOR_PURPLE<<"program"<<COLOR_BLUE<<"] ["<<COLOR_YELLOW<<"arguments"<<COLOR_BLUE<<"]"<<COLOR_DEFAULT<<std::endl;
            cout<<COLOR_BLUE<<"Example: xThermal "<<COLOR_PURPLE<<"thermo"<<COLOR_YELLOW<<" -D0 -VTPX -T100 -P20E5 -H2045E3 -FH2O-NaCl -BIAPWS95 -X0.2"<<COLOR_DEFAULT<<endl;
            cout<<COLOR_BLUE<<"Help:    xThermal ["<<COLOR_PURPLE<<"program"<<COLOR_BLUE<<"] "<<COLOR_YELLOW<<"-h"<<COLOR_BLUE<<COLOR_DEFAULT<<std::endl;
            STATUS("Supported programs are:")
            for(auto & item : m_map_programName2Index)STATUS_color(item.first, COLOR_PURPLE);
        }
        void helpINFO_lutgen()
        {
            STATUS("Help information of "<<COLOR_PURPLE<<"lutgen"<<COLOR_GREEN<<" program");
            STATUS("Generate AMR-LUT of thermodynamic properties and phase regions for specific input fluid.");
            STATUS("[Example]: xThermal lutgen -D2 -VXP -R0.001/0.001/1/1E5/5E5/2200E5 -T600 ");
            int wordWidth=20;
            cout<<"options:"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -h "<<COLOR_BLUE<<"List descriptions of usage and available arguments"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -v "<<COLOR_BLUE<<"Print xThermal version number"<<COLOR_DEFAULT<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -D "<<COLOR_BLUE<<"Dimension: 2, 3. e.g.: -D2 to generate 2D lookup table"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -V "<<COLOR_BLUE<<"Select independent variables according to -D arguments."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Combination of: T, P, X, H. e.g.: -VXT for 2D and -VTPX for 3D"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -T "<<COLOR_BLUE<<"Set fixed temperature value if T is not in -V option."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -P "<<COLOR_BLUE<<"Set fixed pressure value if P is not in -V option. -P316"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -X "<<COLOR_BLUE<<"Set fixed salinity value if X is not in -V option."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -H "<<COLOR_BLUE<<"Set fixed enthalpy value if H is not in -V option."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -R "<<COLOR_BLUE<<"Set range and interval of variables in -V option, must in the save order with -V option."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"e.g.: -R0.001/0.001/1/0.2/1/1000 for -VXT"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -L "<<COLOR_BLUE<<"Set the minimum level for AMR-LUT generation, default value is 4. e.g., -L6."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"     "<<COLOR_DEFAULT<<"Note that the maximum level is automatically calculated from -R options."<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -n "<<COLOR_BLUE<<"If normalize the result in vtk file. Used with -D2 or -D3"<<COLOR_DEFAULT<<std::endl;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -t "<<COLOR_BLUE<<"Set number of thread for parallel computing."<<COLOR_DEFAULT<<std::endl;
        }
        void helpINFO_lutinfo()
        {
            STATUS("Help information of "<<COLOR_PURPLE<<"lutinfo"<<COLOR_GREEN<<" program");
            STATUS("Print AMR-LUT information");
            STATUS("[Usage]: xThermal "<<COLOR_PURPLE<<"lutinfo "<<COLOR_YELLOW<<"myLUT.bin");
        }
        void helpINFO_lut2vtu()
        {
            STATUS("Help information of "<<COLOR_PURPLE<<"lut2vtu"<<COLOR_GREEN<<" program");
            STATUS("Save the AMR-LUT binary file to vtu format for visualization purpose");
            STATUS("[Usage]: xThermal "<<COLOR_PURPLE<<"lut2vtu "<<COLOR_YELLOW<<"myLUT.bin");
            int wordWidth=20;
            cout<<"options:"<<std::endl;;
            cout<<setw(wordWidth)<<setiosflags(ios::left)<<"  -n "<<COLOR_BLUE<<"Normalize [H|P|X|T] as unit length"<<COLOR_DEFAULT<<std::endl;
        }
        void helpINFO()
        {
#ifdef _WIN32
            // set terminal as black(bg)+white(fg) model
      system("color 07"); //see https://www.geeksforgeeks.org/how-to-print-colored-text-in-c/
      GetConsoleScreenBufferInfo(m_hConsole, &csbi);
      m_currentConsoleAttr = csbi.wAttributes;
      int width = (int)(csbi.srWindow.Right-csbi.srWindow.Left+1);
      // int height = (int)(csbi.srWindow.Bottom-csbi.srWindow.Top+1);
      if(width>119)
      {
          StartText_artASCII();
      }else
      {
          StartText();
      }
#else
            struct winsize w;
            ioctl(0, TIOCGWINSZ, &w);
            if(w.ws_col>119)
            {
                StartText_artASCII();
            }else
            {
                StartText();
            }
#endif
            if (m_map_programName2Index.count(m_programName))
            {
                if (m_map_programName2Index[m_programName]==PROGRAM_THERMO)
                {
                    helpINFO_thermo();
                } else if (m_map_programName2Index[m_programName]==PROGRAM_LUTINFO)
                {
                    helpINFO_lutinfo();
                }else if (m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)
                {
                    helpINFO_lutgen();
                }else if (m_map_programName2Index[m_programName]==PROGRAM_LUT2VTU)
                {
                    helpINFO_lut2vtu();
                }else if (m_map_programName2Index[m_programName]==PROGRAM_LOOKUP)
                {
                    helpINFO_lookup();
                }
                else
                {
                    helpINFO_Program();
                }
            }else
            {
                helpINFO_Program();
            }

        }
    private:
        bool GetOptionValue(int opt, char* optarg, double& value);
        // template<typename T>
        // vector<T> linspace(T xmin, T xmax, T dx);

    };
    bool isNum(string str);
}
#endif
