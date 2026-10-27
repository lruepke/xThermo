#ifndef SWEOSBASH_H
#define SWEOSBASH_H
#include "config.h"
#include "getopt_arguments.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <stdexcept>
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
    // Errors of the command-line tool, turned into the exit status by main_cli: 2 for a command line that cannot be
    // interpreted (UsageError), 1 for a failed evaluation (RunError and the library's errors)
    struct UsageError : std::runtime_error { using std::runtime_error::runtime_error; };
    struct RunError : std::runtime_error { using std::runtime_error::runtime_error; };
    #define CLI_USAGE(info) do { std::ostringstream cli_msg_; cli_msg_ << info; throw SWEOSbash::UsageError(cli_msg_.str()); } while (0)
    #define CLI_FAIL(info) do { std::ostringstream cli_msg_; cli_msg_ << info; throw SWEOSbash::RunError(cli_msg_.str()); } while (0)
    #define CLI_WARN(info) do { std::cerr << "Warning: " << info << std::endl; } while (0)

    int main_cli(int argc, char** argv);   // runs the tool and returns its exit status
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
        // Help texts. Every example runs as printed (tests/review/cli/test_cli_help.py); the manual is docs/cli.md.
        static void helpLine(const std::string& option, const std::string& text)
        {
            cout << "  " << setw(6) << setiosflags(ios::left) << option << text << resetiosflags(ios::left) << std::endl;
        }
        static void helpManual()
        {
            cout << "Manual: https://lruepke.github.io/xThermo/cli/" << std::endl;
        }
        static void helpINFO_thermo()
        {
            cout << "xThermal thermo: evaluate states of H2O-NaCl, H2O or NaCl from T-P-X or H-P-X input" << std::endl;
            cout << "Example: xThermal thermo -D0 -VPTX -T300 -P300E5 -X0.1" << std::endl;
            cout << "Example: xThermal thermo -D1 -VT -P250E5 -X0.032 -R1/1/1000 -Oprofile.csv" << std::endl;
            cout << "Options (a value follows the letter directly or after a space: -P300E5 or -P 300E5):" << std::endl;
            helpLine("-D", "0: single points, or a file of points (-G); 1: profile; 2: grid; 3: 3-D grid");
            helpLine("-V", "the varying variables, in order: -D0 PTX or PHX; -D1 T, P, X or H; -D2 PT, PX, TX, PH or HX; -D3 PTX or PHX");
            helpLine("-T", "fixed temperature [deg.C]");
            helpLine("-P", "fixed pressure [Pa], e.g. -P300E5 for 300 bar");
            helpLine("-X", "fixed salinity [kg NaCl/kg fluid, 0-1], e.g. -X0.032 for seawater");
            helpLine("-H", "fixed specific enthalpy [J/kg]; not together with -T");
            helpLine("-R", "min/step/max of each variable of -V, in the order of -V (temperature in deg.C), e.g. -R1/1/1000 for -VT;");
            helpLine("", "both ends are included, the step is adjusted to fit");
            helpLine("-G", "file of points for -D0: three numbers per line in the order of -V, temperature in K");
            helpLine("-O", "output file: CSV for -G and -D1 (required for -D1), VTK (.vtk) for -D2 and -D3");
            helpLine("-F", "fluid: H2O-NaCl (default), H2O, NaCl");
            helpLine("-B", "water model: IAPS84 (default), IAPWS95");
            helpLine("-t", "number of threads; only in a library built with OpenMP (XTHERMO_OPENMP=ON)");
            helpLine("-p", "show a progress bar");
            helpLine("-h", "this help");
            helpLine("-v", "version");
            cout << "Units: temperature in deg.C on the command line, in K in files and output; pressure Pa; salinity kg/kg;" << std::endl;
            cout << "       enthalpy J/kg; density kg/m3; viscosity Pa s. Saturations: volume fractions for H2O-NaCl," << std::endl;
            cout << "       the vapour mass fraction for H2O." << std::endl;
            cout << "Exit status: 0 success, 1 the evaluation failed, 2 the command line cannot be interpreted." << std::endl;
            helpManual();
        }
        void helpINFO_lookup()
        {
            cout << "xThermal lookup: evaluate points from a lookup table built by lutgen" << std::endl;
            cout << "Usage:   xThermal lookup <table.bin> [options]" << std::endl;
            cout << "Example: xThermal lookup lut_constP_XT_6.bin -D0 -VTPX -T300 -P250E5 -X0.1" << std::endl;
            cout << "Options: as for thermo (xThermal thermo -h): -D, -V, -T or -H, -P, -X, -R, -O; the variables and ranges must lie" << std::endl;
            cout << "         in the table. -D0 prints one point, -D1 writes a CSV file, -D2 and -D3 a VTK file." << std::endl;
            helpManual();
        }
        void helpINFO_Program()
        {
            cout << "xThermal " << xThermal_VERSION << " (" << xThermal_DATE << "): thermodynamic properties and phase relations of" << std::endl;
            cout << "H2O-NaCl, H2O and NaCl. xThermo, GEOMAR Helmholtz Centre for Ocean Research Kiel." << std::endl;
            cout << "Usage:   xThermal <program> [options]" << std::endl;
            cout << "Help:    xThermal <program> -h" << std::endl;
            cout << "Example: xThermal thermo -D0 -VPTX -T100 -P20E5 -X0.2 -FH2O-NaCl -BIAPWS95" << std::endl;
            cout << "Programs:" << std::endl;
            helpLine("thermo", "  evaluate states from T-P-X or H-P-X input");
            helpLine("lutgen", "  build an adaptive lookup table");
            helpLine("lutinfo", " describe a lookup table");
            helpLine("lookup", "  evaluate points from a lookup table");
            helpLine("lut2vtu", " convert a lookup table to VTU (ParaView)");
            helpManual();
            cout << "Source: https://github.com/lruepke/xThermo" << std::endl;
        }
        void helpINFO_lutgen()
        {
            cout << "xThermal lutgen: build an adaptive lookup table of the properties over a 2-D or 3-D region" << std::endl;
            cout << "Example: xThermal lutgen -D2 -VXT -R0/0.05/1/100/10/500 -P250E5" << std::endl;
            cout << "Options:" << std::endl;
            helpLine("-D", "2 or 3");
            helpLine("-V", "the varying variables, as for thermo (-D2: PT, PX, TX, PH or HX; -D3: PTX or PHX)");
            helpLine("-T", "fixed temperature [deg.C]; -P pressure [Pa]; -X salinity [kg/kg]; -H enthalpy [J/kg]");
            helpLine("-R", "min/step/max of each variable of -V (temperature in deg.C); the step sets the finest refinement level");
            helpLine("-L", "minimum refinement level (default 4)");
            helpLine("-F", "fluid: H2O-NaCl (default), H2O, NaCl;  -B: water model IAPS84 (default), IAPWS95");
            helpLine("-t", "number of threads; only in a library built with OpenMP (XTHERMO_OPENMP=ON)");
            cout << "Output, in the current directory: lut_<type>_<level>.bin, one file per property, and a .vtu file" << std::endl;
            cout << "(e.g. lut_constP_XT_6.bin for the example). -O is not used." << std::endl;
            helpManual();
        }
        void helpINFO_lutinfo()
        {
            cout << "xThermal lutinfo: describe a lookup table" << std::endl;
            cout << "Usage:   xThermal lutinfo <table.bin>" << std::endl;
            cout << "Example: xThermal lutinfo lut_constP_XT_6.bin" << std::endl;
            helpManual();
        }
        void helpINFO_lut2vtu()
        {
            cout << "xThermal lut2vtu: convert a lookup table to a VTU file for ParaView, <table>.vtu, with normalised axes" << std::endl;
            cout << "Usage:   xThermal lut2vtu <table.bin>" << std::endl;
            cout << "Example: xThermal lut2vtu lut_constP_XT_6.bin" << std::endl;
            helpManual();
        }
        void helpINFO()
        {
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
