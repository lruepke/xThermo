#include "SWEOSbash.h"
#include "LookUpTableForest.h"
#include "LookUpTableForestI.H"
#include <cstdio>
#include <cmath>
#include <set>
// The command-line tool reports errors as exceptions (exit status set in main_cli) and writes errors and warnings to
// stderr. The library's ERROR and WARNING macros print to stdout and call exit(0); they are replaced here.
#undef ERROR
#define ERROR(info) CLI_FAIL(info)
#undef WARNING
#define WARNING(info) CLI_WARN(info)
// using namespace SWEOSbash;
namespace SWEOSbash
{
    int main_cli(int argc, char** argv)
    {
        try
        {
            return bash_run(argc, argv) ? 0 : 1;
        }
        catch (const UsageError& e)
        {
            std::cerr << "Error: " << e.what() << std::endl;
            std::cerr << "Help: xThermal -h, or xThermal [program] -h" << std::endl;
            return 2;
        }
        catch (const std::exception& e)   // RunError and the library's errors
        {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }

    // The output file can be written: checked before a calculation, so that it does not fail at the end
    static void check_writable(const std::string& path)
    {
        const bool existed = std::ifstream(path).good();
        {
            std::ofstream f(path, std::ios::app);
            if (!f) CLI_FAIL("cannot write the output file " << path);
        }
        if (!existed) std::remove(path.c_str());
    }

    // A grid (-D2, -D3) is always written as a legacy VTK file: its name ends in .vtk
    static std::string with_vtk_extension(const std::string& path)
    {
        const size_t slash = path.find_last_of("/\\");
        const size_t dot = path.find_last_of('.');
        const bool hasExtension = dot != std::string::npos && (slash == std::string::npos || dot > slash);
        if (hasExtension && path.substr(dot + 1) == "vtk") return path;
        return (hasExtension ? path.substr(0, dot) : path) + ".vtk";
    }

    bool bash_run(int argc, char** argv)
    {
        // helpINFO();
        //parse arguments and check
        cSWEOSarg arg;
        if(!arg.Parse(argc, argv)) return false;
        if(!arg.Validate(argc, argv))
        {
            // arg.helpINFO();
            return false;
        }
        return true;
    }

    bool isNum(string str) {
        stringstream sin(str);
        double d;
        char c;
        if(!(sin >> d))
            return false;
        if (sin >> c)
            return false;
        return true;
    }

    cSWEOSarg::cSWEOSarg()
            :
            m_pEOS(NULL),
            m_haveD(false), m_haveV(false), m_haveP(false)
            ,m_haveT(false), m_havet(false), m_haveX(false), m_haveH(false), m_haveR(false),m_haveO(false)
            ,m_valueD(-1), m_valueO("out"),m_valueV("")
            ,m_haveB(false), m_valueB("IAPS84")
            ,m_haveG(false), m_haveF(false), m_fluidName("H2O-NaCl")
            ,m_valueX(0.0)
            ,m_normalize_vtk(true)
            ,m_min_level(4), m_max_level(5)
            ,m_TorH(LOOKUPTABLE_FOREST::EOS_ENERGY_T)
            ,m_havep(false)
    {
        m_threadNumOMP = 1;
// #ifndef USE_OMP
//         m_threadNumOMP = 1;
// #else
//         m_threadNumOMP = 1;//omp_get_max_threads();
// #endif
        for(auto & i : m_valueR)
            for(double & j : i)j=0;

        m_map_programName2Index["thermo"] = PROGRAM_THERMO;
        m_map_programName2Index["lutgen"] = PROGRAM_LUTGEN;
        m_map_programName2Index["lutinfo"] = PROGRAM_LUTINFO;
        m_map_programName2Index["lut2vtu"] = PROGRAM_LUT2VTU;
        m_map_programName2Index["lut2pi"] = PROGRAM_LUT2PI;
        m_map_programName2Index["lookup"] = PROGRAM_LOOKUP;

        m_map_fluidName2Ind["H2O"] = FLUID_H2O;
        m_map_fluidName2Ind["H2O-NaCl"] = FLUID_H2ONaCl;
        m_map_fluidName2Ind["NaCl"] = FLUID_NaCl;
        // m_map_fluidName2Ind["Water"] = FLUID_H2O;
        // m_map_fluidName2Ind["water"] = FLUID_H2O;
        // m_map_fluidName2Ind["Salt"] = FLUID_NaCl;
        // m_map_fluidName2Ind["salt"] = FLUID_NaCl;
        // m_map_fluidName2Ind["saltwater"] = FLUID_H2ONaCl;
        // m_map_fluidName2Ind["seawater"] = FLUID_H2ONaCl;
        m_map_backendName2Index["IAPS84"] = BACKEND_H2O_IAPS84;
        m_map_backendName2Index["IAPWS95"] = BACKEND_H2O_IAPWS95;
        m_map_backendName2Index["IAPWS95_CoolProp"] = BACKEND_H2O_IAPWS95_COOLPROP;
    }

    cSWEOSarg::~cSWEOSarg() {

    }

    bool cSWEOSarg::Parse(int argc, char **argv) {
        if(argc<2)
        {
            helpINFO();
            CLI_USAGE("no program given");
        }
        bool needsTableFile = false;
        if(argv[1][0] != '-')
        {
            m_programName = argv[1]; //which program
            optind = 2;// arguments start from the third one
            if (m_map_programName2Index.count(m_programName)>0)
            {
                int program = m_map_programName2Index[m_programName];
                needsTableFile = (program == PROGRAM_LOOKUP || program == PROGRAM_LUTINFO || program == PROGRAM_LUT2VTU || program == PROGRAM_LUT2PI);
            }
            if (needsTableFile && argc>=3 && argv[2][0] != '-')
            {
                m_fileName_LUT = argv[2];
                optind = 3; //arguments start from the fourth one
            }
        } else
        {
            optind = 1;
        }
        int opt;
        const char *optstring = "D:V:P:T:X:H:R:O:G:F:B:t:L:vhnp"; // set argument templete
        int option_index = 0;
        int valid_args=0;
        double doubleOptValue;
        while ((opt = getopt_arguments(argc, argv, optstring)) != -1)
        {
            if(opt=='?') CLI_USAGE("unknown option " << argv[optind-1]);
            if(opt==':') CLI_USAGE("option " << argv[optind-1] << " needs a value");
            valid_args++;
            switch (opt)
            {
                case 'h':
                    helpINFO();
                    exit(0);
                    break;
                case 'v':
                    cout<<"Version: "<<xThermal_VERSION<<endl;
                    exit(0);
                    break;
                case 'n':
                    m_normalize_vtk=true;
                    break;
                case 'D':
                    m_haveD=true;
                    if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
                    if (doubleOptValue != std::floor(doubleOptValue)) CLI_USAGE("-D must be one of 0, 1, 2, 3, not " << optarg);
                    m_valueD=(int)doubleOptValue;
                    break;
                case 't':
                    m_havet=true;
                    if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
#ifdef USE_OMP
                    m_threadNumOMP=(int)doubleOptValue;
                    if(m_threadNumOMP>omp_get_max_threads())m_threadNumOMP=omp_get_max_threads();
                    if(m_threadNumOMP<1)m_threadNumOMP=1;
#endif
                    break;
                case 'V':
                    m_haveV=true;
                    m_valueV=optarg;
                    break;
                case 'P':
                    m_haveP=true;
                    if(!GetOptionValue(opt, optarg, m_valueP))return false;
                    break;
                case 'T':
                    m_haveT=true;
                    if(!GetOptionValue(opt, optarg, m_valueT))return false;
                    break;
                case 'X':
                    m_haveX=true;
                    if(!GetOptionValue(opt, optarg, m_valueX))return false;
                    break;
                case 'H':
                    m_haveH=true;
                    if(!GetOptionValue(opt, optarg, m_valueH))return false;
                    break;
                case 'F':
                    m_haveF = true;
                    m_fluidName = optarg;
                    break;
                case 'G':
                    m_haveG=true;
                    m_valueG=optarg;
                    break;
                case 'R':
                    m_haveR=true;
                    m_valueR_str=xThermal::string_split(optarg,"/");
                    break;
                case 'O':
                    m_haveO=true;
                    m_valueO=optarg;
                    break;
                case 'B': //backend of H2O EOS: [IAPS84|IAPWS95]
                    m_haveB = true;
                    m_valueB = optarg;
                    break;
                case 'L':
                    m_haveL = true;
                    if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
                    m_min_level = max(2, (int)doubleOptValue); //at least set to 2
                    break;
                case 'p':
                    m_havep = true; //show progressbar
                    break;
                default:
                    break;
            }
        }
        if (optind < argc) CLI_USAGE("unexpected argument '" << argv[optind] << "'");
        if(!checkProgram(m_programName))return false;
        if (needsTableFile && m_fileName_LUT.empty())
            CLI_USAGE("program " << m_programName << " needs the file name of a lookup table: xThermal " << m_programName << " table.bin");
        if (!checkFluid(m_fluidName)) return false;

        return true;
    }

    bool cSWEOSarg::GetOptionValue(int opt, char *optarg, double &value) {
        string optarg_str=optarg;
        if(isNum(optarg_str))
        {
            value=atof(optarg);
        }else
        {
            char optCh=opt;
            CLI_USAGE("the value of -"<<optCh<<" is not a number: '"<<optarg_str<<"'");
        }
        return true;
    }

    bool cSWEOSarg::checkFluid(string fluidName) {
#ifdef USE_COOLPROP
        const bool knownBackend = (m_valueB=="IAPS84" || m_valueB=="IAPWS95" || m_valueB=="IAPWS95_CoolProp");
#else
        const bool knownBackend = (m_valueB=="IAPS84" || m_valueB=="IAPWS95");
#endif
        if (!knownBackend) CLI_USAGE("unknown water backend -B" << m_valueB << "; supported: IAPS84, IAPWS95");

        for(auto &item : m_map_fluidName2Ind)
        {
            if(fluidName==item.first)
            {
                if(m_pEOS)
                {
                    delete m_pEOS;
                    m_pEOS = NULL;
                }
                switch (item.second) {
                    case FLUID_H2O:
                        if (m_valueB=="IAPS84")
                        {
                            m_pEOS = new xThermal::PROST::cIAPS84();
                        } else if(m_valueB=="IAPWS95")
                        {
                            m_pEOS = new xThermal::IAPWS95::cIAPWS95();
                        }
#ifdef USE_COOLPROP
                        else if(m_valueB=="IAPWS95_CoolProp")
                        {
                            m_pEOS = new xThermal::COOLPROP::cIAPWS95_CoolProp();
                        }
#endif
                        m_haveX = true;
                        m_valueX = 0;
                    break;
                    case FLUID_NaCl:
                        m_pEOS = new xThermal::NaCl::cNaCl(m_valueB);
                        m_haveX = true;
                        m_valueX = 1.0;
                    break;
                    case FLUID_H2ONaCl:
                        m_pEOS = new xThermal::H2ONaCl::cH2ONaCl(m_valueB);
                    break;
                }
                // set threads num
                m_pEOS->set_num_threads(m_threadNumOMP);
                m_pEOS->showProgressBar(m_havep);
                return true;
            }
        }
        CLI_USAGE("unknown fluid -F" << fluidName << "; supported: H2O-NaCl, H2O, NaCl");
    }

    bool cSWEOSarg::checkProgram(string programName)
    {
        for(auto &item : m_map_programName2Index)
        {
            if(m_programName==item.first)
            {
                return true;
            }
        }
        if (m_programName.empty()) CLI_USAGE("no program given; usage: xThermal [program] [arguments]");
        CLI_USAGE("unknown program '" << m_programName << "'; programs: thermo, lutgen, lutinfo, lut2vtu, lookup");
    }

    bool cSWEOSarg::Validate(int argc, char** argv) {

        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO || m_map_programName2Index[m_programName]==PROGRAM_LUTGEN || m_map_programName2Index[m_programName]==PROGRAM_LOOKUP)
        {
            if(!m_haveD)
            {
                CLI_USAGE("must have -D arguments");
            }
            if(!m_haveV)
            {
                CLI_USAGE("must have -V arguments");
            }
            //check required arguments
            if (m_valueD<0 || m_valueD>3)
            {
                CLI_USAGE("option for -D argument must be one of 0, 1, 2, 3. please check -D parameter");
            }
            for (int i = 0; i < m_valueV.size(); i++)
            {
                if(!(m_valueV[i]=='P' || m_valueV[i]=='T' || m_valueV[i]=='X' || m_valueV[i]=='H'))
                {
                    CLI_USAGE("The option value of -V argument cannot be recognized, the supported variables are combinations of (T, P, X, H), e.g. -VT, -VPX, -VTPX");
                }
            }
            if (m_valueV.size()<1 || m_valueV.size()>3)
            {
                CLI_USAGE("the number of variables cannot exceed three or less than one, it should be one of TPX, THX, T, P, X, H, PT, PX, TX, PH, HX");
            }
            if (m_haveT && m_haveH) CLI_USAGE("give either a temperature (-T) or an enthalpy (-H), not both");
            if (m_haveR)
            {
                // -Rmin/step/max for each variable of -V, in the order of -V
                if (m_valueD == CALCULATION_MODE_SINGLEPOINT) CLI_USAGE("-R is used only with -D1, -D2 and -D3");
                if (m_valueR_str.size() != 3*m_valueV.size())
                    CLI_USAGE("-R needs " << 3*m_valueV.size() << " values (min/step/max for each of " << m_valueV << "), not " << m_valueR_str.size());
                for (size_t i = 0; i < m_valueR_str.size(); i++)
                {
                    if(!isNum(m_valueR_str[i])) CLI_USAGE("value " << i+1 << " of -R is not a number: '" << m_valueR_str[i] << "'");
                    m_valueR[i/3][i%3]=atof(m_valueR_str[i].c_str());
                }
                for (size_t i = 0; i < m_valueV.size(); i++)
                {
                    const double vmin = m_valueR[i][0], step = m_valueR[i][1], vmax = m_valueR[i][2];
                    if (!(vmin < vmax)) CLI_USAGE("-R range of " << m_valueV[i] << ": min (" << vmin << ") must be smaller than max (" << vmax << ")");
                    if (!(step > 0 && step <= vmax - vmin)) CLI_USAGE("-R step of " << m_valueV[i] << " must be > 0 and <= max - min, not " << step);
                }
            }
            if (m_map_programName2Index[m_programName]==PROGRAM_LOOKUP && !std::ifstream(m_fileName_LUT).good())
                CLI_FAIL("cannot open the lookup table " << m_fileName_LUT);
            switch (m_valueD)
            {
                case CALCULATION_MODE_SINGLEPOINT:
                {
                    return Validate_0D();
                }
                    break;
                case CALCULATION_MODE_ONEDIMENSION:
                {
                    return Validate_1D();
                }
                    break;
                case CALCULATION_MODE_TWODIMENSION:
                {
                    return Validate_2D();
                }
                    break;
                case CALCULATION_MODE_THREEDIMENSION:
                {
                    return Validate_3D();
                }
                    break;
                default:
                    break;
            }
        } else if(m_map_programName2Index[m_programName]==PROGRAM_LUTINFO)
        {
            if (!std::ifstream(m_fileName_LUT).good()) CLI_FAIL("cannot open the lookup table " << m_fileName_LUT);
            return Validate_lutinfo(argc, argv);
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LUT2VTU)
        {
            if (!std::ifstream(m_fileName_LUT).good()) CLI_FAIL("cannot open the lookup table " << m_fileName_LUT);
            return Validate_lut2vtu(argc, argv);
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LUT2PI)
        {
            if (!std::ifstream(m_fileName_LUT).good()) CLI_FAIL("cannot open the lookup table " << m_fileName_LUT);
            return Validate_lut2pi(argc, argv);
        }

        return true;
    }

    bool cSWEOSarg::Validate_lutinfo(int argc, char** argv)
    {
        std::string filename_lut = m_fileName_LUT;
        int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(filename_lut);
        if (xThermal::extname_file(filename_lut)!="bin")
        {
            WARNING("The extension name of the input AMR-LUT file is not .bin, please check it since the AMR-LUT file generated by lutgen program is usually saved as .bin file.");
        }
        switch (m_dim_lut)
        {
            case 2:
            {
                LOOKUPTABLE_FOREST::LookUpTableForest_2D lut_2d(filename_lut, NULL, false);
                lut_2d.print_summary();
            }
                break;
            case 3:
            {
                LOOKUPTABLE_FOREST::LookUpTableForest_3D lut_3d(filename_lut, NULL, false);
                lut_3d.print_summary();
            }
                break;
            default:
                ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(filename_lut));
                break;
        }

        return true;
    }

    bool cSWEOSarg::Validate_lut2vtu(int argc, char** argv)
    {
        std::string filename_lut = m_fileName_LUT;
        if (xThermal::extname_file(filename_lut)!="bin")
        {
            WARNING("The extension name of the input AMR-LUT file is not .bin, please check it since the AMR-LUT file generated by lutgen program is usually saved as .bin file.");
        }
        m_pEOS->loadLUT(filename_lut);
        m_pEOS->save_lut_to_vtk(xThermal::filename_without_ext(filename_lut)+".vtu", m_normalize_vtk);
        return true;
    }

    bool cSWEOSarg::Validate_lut2pi(int argc, char** argv)
    {
        std::string filename_lut = m_fileName_LUT;
        if (xThermal::extname_file(filename_lut)!="bin")
        {
            WARNING("The extension name of the input AMR-LUT file is not .bin, please check it since the AMR-LUT file generated by lutgen program is usually saved as .bin file.");
        }
        m_pEOS->loadLUT(filename_lut);
        if(m_pEOS->get_pLUT_lookup())
        {
            if(m_pEOS->get_dim_lut()==2)((LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup())->write_point_index(
                        xThermal::filename_without_ext(filename_lut));
            else ((LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup())->write_point_index(
                        xThermal::filename_without_ext(filename_lut));
        }else
        {
            ERROR("m_pEOS->get_pLUT_lookup() is NULL in bool cSWEOSarg::Validate_lut2pi(int argc, char** argv) function");
        }
        return true;
    }

    bool cSWEOSarg::Validate_0D() {
        if (m_map_programName2Index[m_programName]!=PROGRAM_THERMO && m_map_programName2Index[m_programName]!=PROGRAM_LOOKUP)
        {
            CLI_USAGE("-D0 (single points) is available in the programs thermo and lookup, not in " << m_programName);
        }
        if(m_valueV.size()!=3)
        {
            CLI_USAGE("if -D set as -D0, the -V support [P,T,X] or [P,H,X]: ");
        }
        //single point calculation: PTX
        if(m_valueV=="PTX" || m_valueV=="PXT" || m_valueV=="TPX" || m_valueV=="TXP" || m_valueV=="XPT" || m_valueV=="XTP")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
            if(!m_haveG && !(m_haveT && m_haveP && m_haveX))
                CLI_USAGE("single point " << m_valueV << " needs -T, -P and -X (missing:" << (m_haveT ? "" : " -T") << (m_haveP ? "" : " -P") << (m_haveX ? "" : " -X") << "), or a point file -G");
            // determin single point calculation, multi-points calculation or exit
            if(m_haveT && m_haveP && m_haveX)//single point calculation
            {
                if(m_haveG)CLI_WARN("have -T, -P, -X option values, ignore -G argument");
                //check range
                if(!CheckRange_P(m_valueP))return false;
                if(!CheckRange_T(m_valueT))return false;
                if(!CheckRange_X(m_valueX))return false;
                calculateSinglePoint_PTX(m_valueP, m_valueT+Kelvin, m_valueX);
            }else if(m_haveG)//not T, P, X value, but have intput file name
            {
                STATUS("You specify a input file for multi-points calculation\n"<<COLOR_PURPLE<<"Please make sure your input file with 3 columns in order of "<<m_valueV);
                std::string outputfile = m_valueO;
                if (!m_haveO)
                {
                    outputfile +="_"+m_valueV+".csv";
                    WARNING("The out put file is not specified by -O, I set it to a default filename: "+outputfile);
                }
                check_writable(outputfile);
                calculateMultiPoints_PTX_PHX(m_valueV,m_valueG, outputfile,"T");
            }
        }//single point calculation: PHX
        else if(m_valueV=="PHX" || m_valueV=="PXH" || m_valueV=="HPX" || m_valueV=="HXP" || m_valueV=="XPH" || m_valueV == "XHP")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_H;
            if(!m_haveG && !(m_haveH && m_haveP && m_haveX))
                CLI_USAGE("single point " << m_valueV << " needs -H, -P and -X (missing:" << (m_haveH ? "" : " -H") << (m_haveP ? "" : " -P") << (m_haveX ? "" : " -X") << "), or a point file -G");

            // determin single point calculation, multi-points calculation or exit
            if(m_haveH && m_haveP && m_haveX)//single point calculation
            {
                //check range
                if(!CheckRange_P(m_valueP))return false;
                if(!CheckRange_X(m_valueX))return false;
                calculateSinglePoint_PHX(m_valueP, m_valueH, m_valueX);
            }else if(m_haveG)//not T, P, X value, but have intput file name
            {
                STATUS("You specify a input file for multi-points calculation\n"<<COLOR_PURPLE<<"\tPlease make sure your input file with 3 columns in order of "<<COLOR_RED<<m_valueV);
                std::string outputfile = m_valueO;
                if (!m_haveO)
                {
                    outputfile +="_"+m_valueV+".csv";
                    WARNING("The out put file is not specified by -O, I set it to a default filename: "+outputfile);
                }
                check_writable(outputfile);
                calculateMultiPoints_PTX_PHX(m_valueV,m_valueG, outputfile,"H");
            }
        }
        else
        {
            CLI_USAGE("unknown -V parameter: "<<m_valueV);
        }
        return true;
    }

    bool cSWEOSarg::Validate_1D() {
        if (m_map_programName2Index[m_programName]!=PROGRAM_THERMO && m_map_programName2Index[m_programName]!=PROGRAM_LOOKUP)
        {
            CLI_USAGE("-D1 (profiles) is available in the programs thermo and lookup, not in " << m_programName);
        }
        if(m_valueV.size()!=1)
        {
            CLI_USAGE("if -D set as -D1, the -V option must be one of -VT, -VP, -VX, -VH");
        }
        if (!m_haveR)
        {
            CLI_USAGE("You set -D1 and -V"<<m_valueV<<", then you must set -R for range of "<<m_valueV<<" in format of -Rmin/delta/max(e.g. -R0/1/100)");
        }
        switch (m_valueV[0])
        {
            case 'T'://change T and fixed P, X
            {
                m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
                if(!(m_haveP && CheckRange_P(m_valueP)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument");
                }
                if(!(m_haveX && CheckRange_X(m_valueX)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument");
                }
                if(!m_haveO)
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument");
                }
                check_writable(m_valueO);
                double rangeT[2]={m_valueR[0][0], m_valueR[0][2]}; // [Tmin, Tmax] deg.C
                if(!CheckRanges_T(rangeT)) return false;
                //calculate
                double Trange[3]={m_valueR[0][0]+Kelvin, m_valueR[0][1], m_valueR[0][2]+Kelvin}; // [Tmin, dT, Tmax] K
                vector<double> arrT= m_pEOS->linspace(Trange);
                vector<double> arrP = m_pEOS->linspace(m_valueP, m_valueP, arrT.size());
                vector<double> arrX = m_pEOS->linspace(m_valueX, m_valueX, arrT.size());

                if (m_map_programName2Index[m_programName] == PROGRAM_THERMO)
                {
                    xThermal::ThermodynamicPropertiesVector propsVector = m_pEOS->UpdateState_TPX(arrT, arrP, arrX);
                    propsVector.write(m_valueO);
                }else if (m_map_programName2Index[m_programName] == PROGRAM_LOOKUP)
                {
                    std::vector<std::vector<double>> propsVector;
                    lookup_TorHPX(arrT, arrP, arrX,propsVector, m_valueO);
                }else
                {
                    xThermal::NotImplementedError("Oops! Unexpected case occurs, please check the source code.");
                }
            }
                break;
            case 'P'://change P and fixed T,X or fixed H,X
            {
                //check calculation in TPX space or HPX space
                if(!CheckEOS_Energy())return false;

                if(!(m_haveX && CheckRange_X(m_valueX)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument");
                }
                if(!m_haveO)
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument");
                }
                check_writable(m_valueO);
                double rangeP[2]={m_valueR[0][0], m_valueR[0][2]};
                if(!CheckRanges_P(rangeP)) return false;
                //calculate
                vector<double> arrP= m_pEOS->linspace(m_valueR[0]);
                vector<double> arrTorH = ( m_haveT? m_pEOS->linspace(m_valueT + Kelvin, m_valueT + Kelvin, arrP.size()) : m_pEOS->linspace(m_valueH, m_valueH, arrP.size()));
                vector<double> arrX = m_pEOS->linspace(m_valueX, m_valueX, arrP.size());
                if (m_map_programName2Index[m_programName] == PROGRAM_THERMO)
                {
                    xThermal::ThermodynamicPropertiesVector propsVector = (m_haveT ? m_pEOS->UpdateState_TPX(arrTorH, arrP, arrX) : m_pEOS->UpdateState_HPX(arrTorH, arrP, arrX));
                    propsVector.write(m_valueO);
                }else if (m_map_programName2Index[m_programName] == PROGRAM_LOOKUP)
                {
                    std::vector<std::vector<double>> propsVector;
                    lookup_TorHPX(arrTorH, arrP, arrX,propsVector, m_valueO);
                }else
                {
                    xThermal::NotImplementedError("Oops! Unexpected case occurs, please check the source code.");
                }

            }
                break;
            case 'X'://change X and fixed P, T or P,H
            {
                //check calculation in TPX space or HPX space
                if(!CheckEOS_Energy())return false;
                if(!(m_haveP && CheckRange_P(m_valueP)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument");
                }
                if(!m_haveO)
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument");
                }
                check_writable(m_valueO);
                double rangeX[2]={m_valueR[0][0], m_valueR[0][2]};
                if(!CheckRanges_X(rangeX, 0, 1)) return false;
                //calculate
                vector<double> arrX= m_pEOS->linspace(m_valueR[0]);
                vector<double> arrTorH = ( m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T? m_pEOS->linspace(m_valueT + Kelvin, m_valueT + Kelvin, arrX.size()) : m_pEOS->linspace(m_valueH, m_valueH, arrX.size()));
                vector<double> arrP = m_pEOS->linspace(m_valueP, m_valueP, arrX.size());
                if (m_map_programName2Index[m_programName] == PROGRAM_THERMO)
                {
                    xThermal::ThermodynamicPropertiesVector propsVector = (m_haveT ? m_pEOS->UpdateState_TPX(arrTorH, arrP, arrX) : m_pEOS->UpdateState_HPX(arrTorH, arrP, arrX));
                    propsVector.write(m_valueO);
                }else if (m_map_programName2Index[m_programName] == PROGRAM_LOOKUP)
                {
                    std::vector<std::vector<double>> propsVector;
                    lookup_TorHPX(arrTorH, arrP, arrX,propsVector, m_valueO);
                }else
                {
                    xThermal::NotImplementedError("Oops! Unexpected case occurs, please check the source code.");
                }
            }
                break;
            case 'H'://change H and fixed P, X
            {
                m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_H;
                if(!(m_haveP && CheckRange_P(m_valueP)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument");
                }
                if(!(m_haveX && CheckRange_X(m_valueX)))
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument");
                }
                if(!m_haveO)
                {
                    CLI_USAGE("Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument");
                }
                check_writable(m_valueO);
                // if(!checkranges)
                double rangeH[2]={m_valueR[0][0], m_valueR[0][2]};
                //calculate
                vector<double> arrH= m_pEOS->linspace(m_valueR[0]);
                vector<double> arrP = m_pEOS->linspace(m_valueP, m_valueP, arrH.size());
                vector<double> arrX = m_pEOS->linspace(m_valueX, m_valueX, arrH.size());

                if (m_map_programName2Index[m_programName] == PROGRAM_THERMO)
                {
                    xThermal::ThermodynamicPropertiesVector propsVector = m_pEOS->UpdateState_HPX(arrH, arrP, arrX);
                    propsVector.write(m_valueO);
                }else if (m_map_programName2Index[m_programName] == PROGRAM_LOOKUP)
                {
                    std::vector<std::vector<double>> propsVector;
                    lookup_TorHPX(arrH, arrP, arrX,propsVector, m_valueO);
                }else
                {
                    xThermal::NotImplementedError("Oops! Unexpected case occurs, please check the source code.");
                }
            }
                break;
            default:
                break;
        }

        return true;
    }

    bool cSWEOSarg::Validate_2D() {
        if(m_valueV.size()!=2)
        {
            CLI_USAGE("if -D set as -D2, the -V option must be one of -VPT, -VPX, -VTX, -VPH, -VHX");
        }
        if (!m_haveR)
        {
            CLI_USAGE("You set -D2 and -V"<<m_valueV<<", then you must set -R for range of " << "each variable of -V, as -Rmin/step/max/...");
        }
        if((!m_haveO || m_valueO.empty()) && (m_map_programName2Index[m_programName]!=PROGRAM_LUTGEN))
        {
            CLI_WARN("You forget to set output file name through -O argument, but doesn't matter, it is reseted as "
                <<m_valueV<<".vtk");
            m_valueO=m_valueV+".vtk";
        }
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO && with_vtk_extension(m_valueO) != m_valueO)
        {
            CLI_WARN("grids are written as VTK files: the output file is " << with_vtk_extension(m_valueO) << ", not " << m_valueO);
            m_valueO = with_vtk_extension(m_valueO);
        }
        if (m_map_programName2Index[m_programName]!=PROGRAM_LUTGEN) check_writable(m_valueO);
        // calculate
        vector<double> arrTorH, arrP, arrX;
        const bool cal_meshGrid = true;
        if(m_valueV=="TP" || m_valueV=="PT")// fixed X
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
            if(!(m_haveX && CheckRange_X(m_valueX)))
            {
                CLI_USAGE("Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument");
            }
            int ind_T=0, ind_P=1;
            if(m_valueV=="PT")
            {
                ind_P=0; ind_T=1;
            }
            double rangeT[2]={m_valueR[ind_T][0], m_valueR[ind_T][2]};
            double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
            if(!CheckRanges_T(rangeT)) return false;
            if(!CheckRanges_P(rangeP)) return false;
            //get parameters
            // vector<double> TT, PP;
            // // change unit of T from deg.C to K
            // double Trange_K[3] = {m_valueR[ind_T][0] + Kelvin, m_valueR[ind_T][1], m_valueR[ind_T][2] + Kelvin};
            // std::vector<size_t> nxny= m_pEOS->meshgrid(Trange_K, m_valueR[ind_P], TT, PP);
            // vector<double> XX = m_pEOS->linspace(m_valueX, m_valueX, TT.size());
            // m_pEOS->set_num_threads(8);
            // xThermal::ThermodynamicPropertiesVector propsVector = m_pEOS->UpdateState_TPX(TT,PP,XX);
            // propsVector.write_vtk(m_valueO, m_pEOS->linspace(Trange_K), m_pEOS->linspace(m_valueR[ind_P]), {m_valueX}, "T(K)", "P(Pa)", "X(kg/kg)");

            //change unit of T from deg.C to K
            double Trange_K[3] = {m_valueR[ind_T][0] + Kelvin, m_valueR[ind_T][1], m_valueR[ind_T][2] + Kelvin};
            // calculate max level according to min/delta/max input
            int level_P = ceil(log10((m_valueR[ind_P][2] - m_valueR[ind_P][0])/m_valueR[ind_P][1])/log10(2.0));
            int level_T = ceil(log10((m_valueR[ind_T][2] - m_valueR[ind_T][0])/m_valueR[ind_T][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along P axis: "<<level_P<<", max level along T axis: "<<level_T);
            m_max_level = max(m_min_level, max(m_max_level, max(level_P, level_T))); //maximum level of lutgen
            arrTorH = m_pEOS->linspace(Trange_K);
            arrP = m_pEOS->linspace(m_valueR[ind_P]);
            arrX.push_back(m_valueX);
        }else if(m_valueV=="PX" || m_valueV=="XP")
        {
            //check calculation in TPX space or HPX space
            if(!CheckEOS_Energy())return false;

            int ind_X=0, ind_P=1;
            if(m_valueV=="PX")
            {
                ind_P=0; ind_X=1;
            }
            double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
            double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
            if(!CheckRanges_X(rangeX)) return false;
            if(!CheckRanges_P(rangeP)) return false;
            //get parameters
            // calculate max level according to min/delta/max input
            int level_P = ceil(log10((m_valueR[ind_P][2] - m_valueR[ind_P][0])/m_valueR[ind_P][1])/log10(2.0));
            int level_X = ceil(log10((m_valueR[ind_X][2] - m_valueR[ind_X][0])/m_valueR[ind_X][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along P axis: "<<level_P<<", max level along X axis: "<<level_X);
            m_max_level = max(m_min_level, max(m_max_level, max(level_P, level_X))); //maximum level of lutgen
            //change unit of T from deg.C to K
            arrTorH.push_back(m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? m_valueT + Kelvin : m_valueH);
            arrP = m_pEOS->linspace(m_valueR[ind_P]);
            arrX = m_pEOS->linspace(m_valueR[ind_X]);
        }else if(m_valueV=="TX" || m_valueV=="XT")
        {
            if(!(m_haveP && CheckRange_P(m_valueP)))
            {
                CLI_USAGE("Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument");
            }
            int ind_X=0, ind_T=1;
            if(m_valueV=="TX")
            {
                ind_T=0; ind_X=1;
            }
            double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
            double rangeT[2]={m_valueR[ind_T][0], m_valueR[ind_T][2]};
            if(!CheckRanges_X(rangeX)) return false;
            if(!CheckRanges_T(rangeT)) return false;
            //get parameters
            //change unit of T from deg.C to K
            double Trange_K[3] = {m_valueR[ind_T][0] + Kelvin, m_valueR[ind_T][1], m_valueR[ind_T][2] + Kelvin};
            // calculate max level according to min/delta/max input
            int level_T = ceil(log10((m_valueR[ind_T][2] - m_valueR[ind_T][0])/m_valueR[ind_T][1])/log10(2.0));
            int level_X = ceil(log10((m_valueR[ind_X][2] - m_valueR[ind_X][0])/m_valueR[ind_X][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along T axis: "<<level_T<<", max level along X axis: "<<level_X);
            m_max_level = max(m_min_level, max(m_max_level, max(level_T, level_X))); //maximum level of lutgen
            arrTorH = m_pEOS->linspace(Trange_K);
            arrP.push_back(m_valueP);
            arrX = m_pEOS->linspace(m_valueR[ind_X]);
        }else if(m_valueV=="PH" || m_valueV=="HP")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_H;

            if(!(m_haveX && CheckRange_X(m_valueX)))
            {
                CLI_USAGE("Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument");
            }
            int ind_H=0, ind_P=1;
            if(m_valueV=="PH")
            {
                ind_P=0; ind_H=1;
            }
            double rangeH[2]={m_valueR[ind_H][0], m_valueR[ind_H][2]};
            double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
            if(!CheckRanges_P(rangeP)) return false;
            //get parameters
            // calculate max level according to min/delta/max input
            int level_P = ceil(log10((m_valueR[ind_P][2] - m_valueR[ind_P][0])/m_valueR[ind_P][1])/log10(2.0));
            int level_H = ceil(log10((m_valueR[ind_H][2] - m_valueR[ind_H][0])/m_valueR[ind_H][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along P axis: "<<level_P<<", max level along H axis: "<<level_H);
            m_max_level = max(m_min_level, max(m_max_level, max(level_P, level_H))); //maximum level of lutgen
            //change unit of T from deg.C to K
            arrTorH = m_pEOS->linspace(m_valueR[ind_H]);
            arrP  = m_pEOS->linspace(m_valueR[ind_P]);
            arrX.push_back(m_valueX);
        }else if(m_valueV=="HX" || m_valueV=="XH")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_H;

            if(!(m_haveP && CheckRange_P(m_valueP)))
            {
                CLI_USAGE("Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument");
            }
            int ind_H=0, ind_X=1;
            if(m_valueV=="XH")
            {
                ind_X=0; ind_H=1;
            }
            double rangeH[2]={m_valueR[ind_H][0], m_valueR[ind_H][2]};
            double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
            if(!CheckRanges_X(rangeX)) return false;
            //get parameters
            // calculate max level according to min/delta/max input
            int level_X = ceil(log10((m_valueR[ind_X][2] - m_valueR[ind_X][0])/m_valueR[ind_X][1])/log10(2.0));
            int level_H = ceil(log10((m_valueR[ind_H][2] - m_valueR[ind_H][0])/m_valueR[ind_H][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along X axis: "<<level_X<<", max level along H axis: "<<level_H);
            m_max_level = max(m_min_level, max(m_max_level, max(level_X, level_H))); //maximum level of lutgen
            //change unit of T from deg.C to K
            arrTorH = m_pEOS->linspace(m_valueR[ind_H]);
            arrP.push_back(m_valueP);
            arrX = m_pEOS->linspace(m_valueR[ind_X]);
        }else
        {
            CLI_USAGE("-V" << m_valueV << " is not available with -D2; use one of PT, PX, TX, PH, HX");
        }

        //calculate
        int UpdateProps_lutGen = Update_prop_CommonlyUsedHPX;
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO)
        {
            xThermal::ThermodynamicPropertiesVector propsVector = (m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H ? m_pEOS->UpdateState_HPX(arrTorH,arrP,arrX, cal_meshGrid) : m_pEOS->UpdateState_TPX(arrTorH,arrP,arrX, cal_meshGrid));
            propsVector.write_vtk(m_valueO, arrX, arrTorH, arrP, "X(kg/kg)", (m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H ? "H(J/kg)" : "T(K)"), "P(Pa)");
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)
        {
            if (arrP.size()==1) //const P
            {
                m_pEOS->createLUT_2D(*std::min_element(arrX.begin(), arrX.end()), *std::max_element(arrX.begin(), arrX.end()),
                                     *std::min_element(arrTorH.begin(), arrTorH.end()),*std::max_element(arrTorH.begin(), arrTorH.end()),
                                     arrP[0], LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, m_TorH, m_min_level, m_max_level, UpdateProps_lutGen);
                m_pEOS->save_lut_to_binary((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constP_XT_" : "lut_constP_XH_")+std::to_string(m_max_level)+".bin");
                m_pEOS->save_lut_to_vtk((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constP_XT_" : "lut_constP_XH_")+std::to_string(m_max_level)+".vtu", m_normalize_vtk);
            }else if (arrTorH.size()==1) //const TorH
            {
                m_pEOS->createLUT_2D(*std::min_element(arrX.begin(), arrX.end()), *std::max_element(arrX.begin(), arrX.end()),
                                     *std::min_element(arrP.begin(), arrP.end()),*std::max_element(arrP.begin(), arrP.end()),
                                     arrTorH[0], LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, m_TorH, m_min_level, m_max_level, UpdateProps_lutGen);
                m_pEOS->save_lut_to_binary((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constT_XP_" : "lut_constH_XP_")+std::to_string(m_max_level)+".bin");
                m_pEOS->save_lut_to_vtk((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constT_XP_" : "lut_constH_XP_")+std::to_string(m_max_level)+".vtu", m_normalize_vtk);
            }else if (arrX.size()==1) //const X
            {
                m_pEOS->createLUT_2D(*std::min_element(arrTorH.begin(), arrTorH.end()), *std::max_element(arrTorH.begin(), arrTorH.end()),
                                     *std::min_element(arrP.begin(), arrP.end()),*std::max_element(arrP.begin(), arrP.end()),
                                     arrX[0], LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, m_TorH, m_min_level, m_max_level, UpdateProps_lutGen);
                m_pEOS->save_lut_to_binary((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constX_TP_" : "lut_constX_HP_")+std::to_string(m_max_level)+".bin");
                m_pEOS->save_lut_to_vtk((m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "lut_constX_TP_" : "lut_constX_HP_")+std::to_string(m_max_level)+".vtu", m_normalize_vtk);
            }
        }else if(m_map_programName2Index[m_programName] == PROGRAM_LOOKUP)
        {
            std::vector<std::vector<double>> props;
            lookup_TorHPX(arrTorH, arrP, arrX,props, m_valueO, cal_meshGrid);
        }

        return true;
    }

    bool cSWEOSarg::Validate_3D() {
        if (m_map_fluidName2Ind[m_fluidName]!=FLUID_H2ONaCl)
        {
            CLI_USAGE("-D3 is available only for the fluid H2O-NaCl, not " << m_fluidName);
        }
        if(m_valueV.size()!=3)
        {
            CLI_USAGE("if -D set as -D3, the -V option must be one of -VPTX, -VPHX");
        }
        if (!m_haveR)
        {
            CLI_USAGE("You set -D3 and -V"<<m_valueV<<", then you must set -R for range of " << "each variable of -V, as -Rmin/step/max/...");
        }
        if((!m_haveO || m_valueO.empty()) && (m_map_programName2Index[m_programName]!=PROGRAM_LUTGEN))
        {
            CLI_WARN("You forget to set output file name through -O argument, but doesn't matter, it is reseted as "
                <<m_valueV<<".vtk");
            m_valueO=m_valueV+".vtk";
        }
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO && with_vtk_extension(m_valueO) != m_valueO)
        {
            CLI_WARN("grids are written as VTK files: the output file is " << with_vtk_extension(m_valueO) << ", not " << m_valueO);
            m_valueO = with_vtk_extension(m_valueO);
        }
        if (m_map_programName2Index[m_programName]!=PROGRAM_LUTGEN) check_writable(m_valueO);
        //calculate
        vector<double> arrTorH, arrP, arrX;
        const bool cal_meshGrid = true;
        if(m_valueV=="PTX" || m_valueV=="PXT" || m_valueV=="TPX" || m_valueV=="TXP" || m_valueV=="XPT" || m_valueV=="XTP")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
            int indP=0, indT=1, indX=2;
            if(m_valueV=="PXT")
            {
                indP=0; indX=1; indT=2;
            }else if(m_valueV=="TPX")
            {
                indT=0; indP=1; indX=2;
            }else if(m_valueV=="TXP")
            {
                indT=0; indX=1; indP=2;
            }else if(m_valueV=="XPT")
            {
                indX=0; indP=1; indT=2;
            }else if(m_valueV=="XTP")
            {
                indX=0; indT=1; indP=2;
            }
            double rangeT[2]={m_valueR[indT][0], m_valueR[indT][2]};
            double rangeP[2]={m_valueR[indP][0], m_valueR[indP][2]};
            double rangeX[2]={m_valueR[indX][0], m_valueR[indX][2]};
            if(!CheckRanges_T(rangeT)) return false;
            if(!CheckRanges_P(rangeP)) return false;
            if(!CheckRanges_X(rangeX)) return false;
            //calculate
            //change unit of T from deg.C to K
            double Trange_K[3] = {m_valueR[indT][0] + Kelvin, m_valueR[indT][1], m_valueR[indT][2] + Kelvin};
            // calculate max level according to min/delta/max input
            int level_P = ceil(log10((m_valueR[indP][2] - m_valueR[indP][0])/m_valueR[indP][1])/log10(2.0));
            int level_T = ceil(log10((m_valueR[indT][2] - m_valueR[indT][0])/m_valueR[indT][1])/log10(2.0));
            int level_X = ceil(log10((m_valueR[indX][2] - m_valueR[indX][0])/m_valueR[indX][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along P axis: "<<level_P<<", max level along T axis: "<<level_T<<", max level along X axis: "<<level_X);
            m_max_level = max(m_min_level, max(m_max_level, max(level_P, max(level_T, level_X)))); //maximum level of lutgen
            arrTorH = m_pEOS->linspace(Trange_K);
            arrP = m_pEOS->linspace(m_valueR[indP]);
            arrX = m_pEOS->linspace(m_valueR[indX]);
        }else if(m_valueV=="PHX" || m_valueV=="PXH" || m_valueV=="HPX" || m_valueV=="HXP" || m_valueV=="XPH" || m_valueV=="XHP")
        {
            m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_H;
            int indP=0, indH=1, indX=2;
            if(m_valueV=="PXH")
            {
                indP=0; indX=1; indH=2;
            }else if(m_valueV=="HPX")
            {
                indH=0; indP=1; indX=2;
            }else if(m_valueV=="HXP")
            {
                indH=0; indX=1; indP=2;
            }else if(m_valueV=="XPH")
            {
                indX=0; indP=1; indH=2;
            }else if(m_valueV=="XHP")
            {
                indX=0; indH=1; indP=2;
            }
            double rangeH[2]={m_valueR[indH][0], m_valueR[indH][2]};
            double rangeP[2]={m_valueR[indP][0], m_valueR[indP][2]};
            double rangeX[2]={m_valueR[indX][0], m_valueR[indX][2]};
            double rangePX[4]={m_valueR[indP][0], m_valueR[indP][2], m_valueR[indX][0], m_valueR[indX][2]};
            if(!CheckRanges_P(rangeP)) return false;
            if(!CheckRanges_X(rangeX)) return false;
            //calculate
            // calculate max level according to min/delta/max input
            int level_P = ceil(log10((m_valueR[indP][2] - m_valueR[indP][0])/m_valueR[indP][1])/log10(2.0));
            int level_H = ceil(log10((m_valueR[indH][2] - m_valueR[indH][0])/m_valueR[indH][1])/log10(2.0));
            int level_X = ceil(log10((m_valueR[indX][2] - m_valueR[indX][0])/m_valueR[indX][1])/log10(2.0));
            if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)STATUS("max level along P axis: "<<level_P<<", max level along T axis: "<<level_H<<", max level along X axis: "<<level_X);
            m_max_level = max(m_min_level, max(m_max_level, max(level_P, max(level_H, level_X)))); //maximum level of lutgen
            arrTorH= m_pEOS->linspace(m_valueR[indH]);
            arrP = m_pEOS->linspace(m_valueR[indP]);
            arrX = m_pEOS->linspace(m_valueR[indX]);
        }

        //calculate
        int UpdateProps_lutGen = Update_prop_CommonlyUsedHPX;
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO)
        {
            xThermal::ThermodynamicPropertiesVector propsVector = ( m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H ? m_pEOS->UpdateState_HPX(arrTorH,arrP,arrX, cal_meshGrid) : m_pEOS->UpdateState_TPX(arrTorH,arrP,arrX, cal_meshGrid) );
            propsVector.write_vtk(m_valueO, arrX, arrTorH, arrP, "X(kg/kg)", (m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H ? "H(J/kg)" : "T(K)"), "P(Pa)");
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LUTGEN)
        {
            if (m_max_level>15)
            {
                WARNING("The input max level for the 3D AMR-LUT is greater than 15, which will generate a huge table and need huge amount of memory and disk space, do you want continue ?");
                std::string yes_no;
                std::cout<<"[y|n]?";
                std::cin>>yes_no;
                if(yes_no!="y")
                {
                    STATUS("You do not select [y]es, so stop calculating.");
                    exit(0);
                }
            }

            m_pEOS->createLUT_3D(*std::min_element(arrTorH.begin(), arrTorH.end()), *std::max_element(arrTorH.begin(), arrTorH.end()),
                                 *std::min_element(arrP.begin(), arrP.end()),*std::max_element(arrP.begin(), arrP.end()),
                                 *std::min_element(arrX.begin(), arrX.end()),*std::max_element(arrX.begin(), arrX.end()),
                                 m_TorH, m_min_level, m_max_level, UpdateProps_lutGen);
            m_pEOS->save_lut_to_binary("lut_"+m_valueV+"_"+std::to_string(m_max_level)+".bin");
            m_pEOS->save_lut_to_vtk("lut_"+m_valueV+"_"+std::to_string(m_max_level)+".vtu", m_normalize_vtk);
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LOOKUP)
        {
            std::vector<std::vector<double>> props;
            lookup_TorHPX(arrTorH, arrP, arrX,props, m_valueO, cal_meshGrid);
        }

        return true;
    }

    bool cSWEOSarg::CheckRange_T(double T_C) {
        return CheckRange_TorH(LOOKUPTABLE_FOREST::EOS_ENERGY_T, T_C + Kelvin, m_pEOS->Tmin(), m_pEOS->Tmax(), "");
    }

    bool cSWEOSarg::CheckRange_P(double P0) {
        return CheckRange_P(P0, m_pEOS->pmin(), m_pEOS->pmax());
    }

    bool cSWEOSarg::CheckRange_P(double P0, double Pmin, double Pmax, const std::string& info) {
        if(P0>Pmax || P0<Pmin)
        {
            CLI_FAIL("-P specify pressure ="<<P0/1E5<<"bar = "<<P0<<" Pa, out of range ["<<Pmin<<", "<<Pmax<<"] Pa "<<info);
        }
        return true;
    }

    bool cSWEOSarg::CheckRange_TorH(LOOKUPTABLE_FOREST::EOS_ENERGY eos_energy,double TorH0, double TorHmin, double TorHmax, const std::string& info) {
        switch (eos_energy) {
            case LOOKUPTABLE_FOREST::EOS_ENERGY_T:
            {
                if(TorH0>TorHmax || TorH0<TorHmin)
                {
                    CLI_FAIL("-T specify temperature ="<<TorH0-Kelvin<<"deg.C = "<<TorH0<<" K, out of range ["<<TorHmin<<", "<<TorHmax<<"] K "<<info);
                }
            }
                break;
            case LOOKUPTABLE_FOREST::EOS_ENERGY_H:
            {
                if(TorH0>TorHmax || TorH0<TorHmin)
                {
                    CLI_FAIL("-H specify enthalpy ="<<TorH0/1000<<" kJ/kg = "<<TorH0<<" J/kg, out of range ["<<TorHmin<<", "<<TorHmax<<"] J/kg "<<info);
                }
            }
                break;
        }
        return true;
    }

    bool cSWEOSarg::CheckRange_X(double X0, double XMIN, double XMAX, const std::string& info) {
        if(X0>XMAX || X0<XMIN)
        {
            CLI_FAIL("-X specify salinity ="<<X0<<" wt. % NaCl, out of range ["<<XMIN<<", "<<XMAX<<"] wt. % NaCl "<<info);
        }
        return true;
    }

    xThermal::ThermodynamicProperties cSWEOSarg::calculateSinglePoint_PTX(double P, double T_K, double X, bool isCout) {
        xThermal::ThermodynamicProperties prop;
        prop.fluidName = m_pEOS->name();
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO)
        {
            prop = m_pEOS->UpdateState_TPX(T_K, P, X);
            if(isCout)
            {
                cout<<prop<<endl;
            }
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LOOKUP)
        {
            bool printStatus_lut = false;
            m_pEOS->loadLUT(m_fileName_LUT,printStatus_lut);
            const int dim = m_pEOS->get_dim_lut_lookup();
            int num_props_per_node = 0;
            switch (dim) {
                case 2:
                {
                    auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
                    num_props_per_node = (int)pLUT->m_map_props.size();
                    std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                    // std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                    if (pLUT->m_TorH != LOOKUPTABLE_FOREST::EOS_ENERGY_T) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 2D table in "<<name_space<<" space, but your input -V is "<<m_valueV);
                }
                    break;
                case 3:
                {
                    auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
                    num_props_per_node = (int)pLUT->m_map_props.size();
                    std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                    std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                    if (pLUT->m_TorH != LOOKUPTABLE_FOREST::EOS_ENERGY_T) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 3D table in "<<name_space<<" space, but your input -V is "<<m_valueV);
                }
                    break;
                default: ERROR("Something is wrong in the AMR-LUT, because the dim is neither 2 nor 3, please check the LUT file: "<<m_fileName_LUT);
                    break;
            }
            double props[num_props_per_node];
            lookup_TorHPX(dim, T_K, P, X, props, isCout);
        }
        return prop;
    }

    xThermal::PhaseRegion cSWEOSarg::lookup_TorHPX(const int& dim, const double& TorH, const double& P, const double& X, double* props, bool isCout) {
        xThermal::PhaseRegion phaseRegion;
        // check input parameters
        if (dim==2)
            {
                auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
                double xyz_min_target[dim];
                switch (pLUT->m_const_which_var) {
                    case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                    {
                        CheckRange_TorH(pLUT->m_TorH,TorH, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT);//index of T is zero according to CONST_X_VAR_TorHP
                        CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT);//index of P is zero according to CONST_X_VAR_TorHP
                        auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, TorH, P);
                        phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                        if (isCout)
                        {
                            if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, P = "<<P<<" Pa")
                            else STATUS("Input: H = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, P = "<<P<<" Pa");
                            STATUS("Constant in the 2D LUT: X = "<<pLUT->m_constZ<<" kg/kg");
                        }
                    }
                        break;
                    case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                    {
                        CheckRange_X(X, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT);//index of T is zero according to CONST_X_VAR_TorHP
                        CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT);//index of P is zero according to CONST_X_VAR_TorHP
                        auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, X, P);
                        phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                        if (isCout)
                        {
                            STATUS("Input: P = "<<P<<" Pa, X = "<<X<<" kg/kg");
                            if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Constant in the 2D LUT: T = "<<pLUT->m_constZ-Kelvin<<" deg.C = "<<pLUT->m_constZ<<" K")
                            else STATUS("Constant in the 2D LUT: T = "<<pLUT->m_constZ/1000<<" kJ/kg = "<<pLUT->m_constZ<<" J/kg");
                        }
                    }
                        break;
                    case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                    {
                        CheckRange_X(X, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT);//index of T is zero according to CONST_X_VAR_TorHP
                        CheckRange_TorH(pLUT->m_TorH, TorH, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT);//index of P is zero according to CONST_X_VAR_TorHP
                        auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, X, TorH);
                        phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                        if (isCout)
                        {
                            if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, X = "<<X<<" kg/kg")
                            else STATUS("Input: H = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, X = "<<X<<" kg/kg");
                            STATUS("Constant in the 2D LUT: P = "<<pLUT->m_constZ<<" Pa");
                        }
                    }
                        break;
                    default: ERROR("pLUT->m_const_which_var of 2D AMR-LUT is wrong: "<<pLUT->m_const_which_var<<", it is impossible. Check source code: ThermodynamicProperties cSWEOSarg::calculateSinglePoint_PTX(double P, double T_K, double X, bool isCout)");
                        break;
                }
                if (isCout)
                {
                    int i = 0;
                    STATUS("Properties lookup result");
                    cout<<"    "<<COLOR_BLUE<<"Phase region: "<<COLOR_PURPLE<<xThermal::map_phase2name[phaseRegion]<<COLOR_DEFAULT<<endl;
                    for(auto &item : pLUT->m_map_props)
                    {
                        cout<<"    "<<COLOR_BLUE<<item.second.longName<<COLOR_PURPLE<<": "<<props[i]<<COLOR_DEFAULT<<" "<<item.second.unit<<endl;
                        i++;
                    }
                }
            }else if(dim==3)
            {
                auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
                // double props[pLUT->m_map_props.size()];
                double xyz_min_target[dim];
                CheckRange_TorH(pLUT->m_TorH, TorH, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT);//index of T is zero according to CONST_X_VAR_TorHP
                CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT);//index of P is zero according to CONST_X_VAR_TorHP
                CheckRange_X(X, pLUT->m_xyz_min[2], pLUT->m_xyz_max[2], "in AMR-LUT "+m_fileName_LUT);//index of P is zero according to CONST_X_VAR_TorHP
                auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, TorH, P, X);
                phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                if (isCout)
                {
                    if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, P = "<<P<<" Pa, X = "<<X<<" kg/kg")
                    else STATUS("Input: T = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, P = "<<P<<" Pa, X = "<<X<<" kg/kg");
                    int i = 0;
                    STATUS("Properties lookup result");
                    cout<<"    "<<COLOR_BLUE<<"Phase region: "<<COLOR_PURPLE<<xThermal::map_phase2name[phaseRegion]<<COLOR_DEFAULT<<endl;
                    for(auto &item : pLUT->m_map_props)
                    {
                        cout<<"    "<<COLOR_BLUE<<item.second.longName<<COLOR_PURPLE<<": "<<props[i]<<COLOR_DEFAULT<<" "<<item.second.unit<<endl;
                        i++;
                    }
                }
            }
        return phaseRegion;
    }

    /**
    *
    * @param TorH [SI]
    * @param P [SI]
    * @param X [SI]
    * @param props 2D vector, N rows and M columns. M is property number, the same as size of the reture vector. Note that the last column of each row is the phase region index.
    * @param filename_out_csv
    * @param isMeshGrid
    * @return Property names vector.
    */
    std::vector<xThermal::propInfo> cSWEOSarg::lookup_TorHPX(const vector<double> &TorH, const vector<double> &P, const vector<double> &X,
                                  std::vector<std::vector<double>>& props, std::string filename_out, bool isMeshGrid )
    {
        bool printStatus_lut = false;
        m_pEOS->loadLUT(m_fileName_LUT,printStatus_lut);
        const int dim = m_pEOS->get_dim_lut_lookup();
        int num_props_per_node = 0;
        std::string name_props;
        std::vector<xThermal::propInfo> propsInfo;
        switch (dim) {
            case 2:
            {
                auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
                num_props_per_node = (int)pLUT->m_map_props.size();
                std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                // std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                if (pLUT->m_TorH != m_TorH) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 2D table in "<<name_space<<" space, but your input TPX/HPX space is not match with that of the LUT file.");
                for(auto &item : pLUT->m_map_props)
                {
                    propsInfo.push_back(item.second);
                    name_props += "," + std::string(item.second.shortName)+ "(" + std::string(item.second.unit) + ")";
                }
            }
                break;
            case 3:
            {
                auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
                num_props_per_node = (int)pLUT->m_map_props.size();
                std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                if (pLUT->m_TorH != m_TorH) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 3D table in "<<name_space<<" space, but your input TPX/HPX space is not match with that of the LUT file.");
                for(auto &item : pLUT->m_map_props)
                {
                    propsInfo.push_back(item.second);
                    name_props += "," + std::string(item.second.shortName)+ "(" + std::string(item.second.unit) + ")";
                }
            }
                break;
            default: ERROR("Something is wrong in the AMR-LUT, because the dim is neither 2 nor 3, please check the LUT file: "<<m_fileName_LUT);
                break;
        }
        xThermal::PhaseRegion phaseRegion;
        double props_[num_props_per_node];
        propsInfo.push_back({"Phase","Phase region","-"});
        //lookup and save to file
        if (isMeshGrid)
        {
            size_t nTorH = TorH.size(), nP = P.size(), nX = X.size();
            size_t N = nTorH*nP*nX;
            size_t nXnTorH = nX*nTorH;
            std::vector<double> TorH_vec(N), P_vec(N), X_vec(N);
            props.clear();
            props.resize(N);
            int jj = 0;
            // stateVector.resize(N);
            // ThermodynamicProperties props;
            // size_t  ind = 0;
            MultiProgressBar multiBar(nTorH*nP);
// #ifdef USE_OMP
//             if(m_pEOS->get_num_threads()>1)STATUS("Parallel computing, threads number: "<<m_pEOS->get_num_threads()<<"\n");
// #pragma omp parallel for private(props_, m_pEOS) shared(props, TorH, P, X)
// #endif
            for (int i = 0; i < nP; ++i) {
                for (int j = 0; j < nTorH; ++j) {
                    for (int k = 0; k < nX; ++k) {
                        size_t  ind = k + j*nX + i*(nXnTorH);
                        TorH_vec[ind] = TorH[j];    P_vec[ind] = P[i];  X_vec[ind] = X[k];
                        phaseRegion = lookup_TorHPX(dim, TorH[j], P[i], X[k], props_);
                        props[ind].resize(num_props_per_node + 1);//the last one is phase index
                        for (jj = 0; jj < num_props_per_node; ++jj) {
                            props[ind][jj] = props_[jj];
                        }
                        props[ind][jj] = phaseRegion;
                    }
// #ifdef USE_OMP
// #pragma omp critical
// #endif
                    multiBar.Update();
                }
            }
            // save to file
            if (!filename_out.empty())
            {
                if (xThermal::extname_file(filename_out)!="vtk")
                {
                    filename_out = xThermal::filename_without_ext(filename_out)+".vtk";
                    WARNING("The output file name is not .vtk file, I will change it to "<<COLOR_BLUE<<filename_out);
                }
                m_pEOS->writeMeshGrid2VTK(filename_out, X, "X(kg/kg)", TorH, (m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H ? "H(J/kg)" : "T(K)"), P, "P(Pa)", props, propsInfo);
            }

        } else
        {
            //clear props vector before filling it.
            props.clear();
            props.resize(TorH.size());
            int j =0;
            for (int i = 0; i < TorH.size(); ++i) {
                phaseRegion = lookup_TorHPX(dim, TorH[i], P[i], X[i], props_);
                props[i].resize(num_props_per_node + 1);//the last one is phase index
                for (j = 0; j < num_props_per_node; ++j) {
                    props[i][j] = props_[j];
                }
                props[i][j] = phaseRegion;
            }
            //save to file
            if (!filename_out.empty())
            {
                if (xThermal::extname_file(filename_out)!="csv")
                {
                    filename_out = xThermal::filename_without_ext(filename_out)+".csv";
                    WARNING("The output file name is not .csv file, I will change it to "<<COLOR_BLUE<<filename_out);
                }
                ofstream fpout(filename_out);
                if (!fpout) ERROR("Open file failed in cSWEOSarg::lookup_TorHPX(const vector<double> &TorH, const vector<double> &P, const vector<double> &X, ...): "<<filename_out);
                fpout<<(m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T(K)" : "H(J/kg)")<<", P(Pa), X(kg/kg)"<<name_props<<", Phase index"<<endl;
                size_t num_props_total = props[0].size();
                for (int i = 0; i < TorH.size(); ++i) {
                    fpout<<TorH[i]<<", "<<P[i]<<", "<<X[i];
                    for (j = 0; j < num_props_total; ++j) {
                        fpout<<", "<<props[i][j];
                    }
                    fpout<<endl;
                }
                fpout.close();
            }
        }
        return propsInfo;
    }

    xThermal::ThermodynamicProperties cSWEOSarg::calculateSinglePoint_PHX(double P, double H, double X, bool isCout) {
        xThermal::ThermodynamicProperties prop;
        prop.fluidName = m_pEOS->name();
        if (m_map_programName2Index[m_programName]==PROGRAM_THERMO)
        {
            prop = m_pEOS->UpdateState_HPX(H, P, X);
            if(isCout)
            {
                cout<<prop<<endl;
            }
        }else if(m_map_programName2Index[m_programName]==PROGRAM_LOOKUP)
        {
            bool printStatus_lut = false;
            m_pEOS->loadLUT(m_fileName_LUT,printStatus_lut);
            const int dim = m_pEOS->get_dim_lut_lookup();
            int num_props_per_node = 0;
            switch (dim) {
                case 2:
                {
                    auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
                    num_props_per_node = (int)pLUT->m_map_props.size();
                    std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                    // std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                    if (pLUT->m_TorH != LOOKUPTABLE_FOREST::EOS_ENERGY_H) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 2D table in "<<name_space<<" space, but your input -V is "<<m_valueV);
                }
                    break;
                case 3:
                {
                    auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
                    num_props_per_node = (int)pLUT->m_map_props.size();
                    std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
                    std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
                    if (pLUT->m_TorH != LOOKUPTABLE_FOREST::EOS_ENERGY_H) ERROR("The input AMR-LUT: "<<m_fileName_LUT<<" is a 3D table in "<<name_space<<" space, but your input -V is "<<m_valueV);
                }
                    break;
                default: ERROR("Something is wrong in the AMR-LUT, because the dim is neither 2 nor 3, please check the LUT file: "<<m_fileName_LUT);
                    break;
            }
            double props[num_props_per_node];
            lookup_TorHPX(dim, H, P, X, props, isCout);

        }
        return prop;
    }

        xThermal::ThermodynamicPropertiesVector
    cSWEOSarg::calculateMultiPoints_PTX_PHX(const string& valueV, const string& filePTX, string outFile, const string& isT_H) {
        ifstream fin(filePTX);
        if(!fin) CLI_FAIL("cannot open the point file -G" << filePTX);
        // column order of the file: the order of -V (T in K, P in Pa, X in kg/kg, H in J/kg)
        const std::string TorH = isT_H;
        std::set<std::string> orders = {"P"+TorH+"X", "PX"+TorH, TorH+"PX", TorH+"XP", "XP"+TorH, "X"+TorH+"P"};
        if (orders.count(valueV) == 0) CLI_USAGE("-V" << valueV << " is not a column order of P, " << TorH << " and X");
        vector<double>P, T_H, X;
        std::string line;
        size_t lineNumber = 0;
        while (std::getline(fin, line))
        {
            lineNumber++;
            if (line.find_first_not_of(" \t\r") == std::string::npos) continue;   // blank line
            std::istringstream fields(line);
            double v[3];
            std::string extra;
            if (!(fields >> v[0] >> v[1] >> v[2]) || (fields >> extra))
                CLI_FAIL("point file " << filePTX << ", line " << lineNumber << ": expected three numbers in the order " << valueV << ", got '" << line << "'");
            for (int k = 0; k < 3; k++)
            {
                if (valueV[k] == 'P') P.push_back(v[k]);
                else if (valueV[k] == 'X') X.push_back(v[k]);
                else T_H.push_back(v[k]);
            }
        }
        if (P.empty()) CLI_FAIL("point file " << filePTX << " contains no points");
        fin.close();
        //calculate
        xThermal::ThermodynamicPropertiesVector propsVector;
        // MultiProgressBar multibar(P.size(),COLOR_BAR_BLUE);
        if(isT_H=="T")
        {
            propsVector = m_pEOS->UpdateState_TPX(T_H,P,X);
        }else if(isT_H=="H")
        {
            propsVector = m_pEOS->UpdateState_HPX(T_H,P,X);
        }
        //output
        propsVector.write(outFile);
        return propsVector;
    }

    bool cSWEOSarg::CheckRanges_T(double *Trange) {
        if(Trange[0]>m_pEOS->Tmax()-Kelvin || Trange[0]<m_pEOS->Tmin()-Kelvin)
        {
            CLI_FAIL("The minimum value of temperature specified by -R argument Tmin ="<<Trange[0]<<"deg.C = "<<Trange[0]+Kelvin<<" K, out of range ["<<m_pEOS->Tmin()<<", "<<m_pEOS->Tmax()<<"] K");
        }
        if(Trange[1]>m_pEOS->Tmax()-Kelvin || Trange[1]<m_pEOS->Tmin()-Kelvin)
        {
            CLI_FAIL("The maximum value of temperature specified by -R argument Tmax ="<<Trange[1]<<"deg.C = "<<Trange[1]+Kelvin<<" K, out of range ["<<m_pEOS->Tmin()<<", "<<m_pEOS->Tmax()<<"] K");
        }
        return true;
    }

    bool cSWEOSarg::CheckRanges_P(double *Prange) {
        if(Prange[0]>m_pEOS->pmax() || Prange[0]<m_pEOS->pmin())
        {
            CLI_FAIL("The minimum value of pressure specified by -R argument Pmin ="<<Prange[0]<<"bar = "<<Prange[0]*1e5<<" Pa, out of range ["<<m_pEOS->pmin()<<", "<<m_pEOS->pmax()<<"] Pa");
        }
        if(Prange[1]>m_pEOS->pmax() || Prange[1]<m_pEOS->pmin())
        {
            CLI_FAIL("The minimum value of pressure specified by -R argument Pmax ="<<Prange[1]<<"bar = "<<Prange[1]*1e5<<" Pa, out of range ["<<m_pEOS->pmin()<<", "<<m_pEOS->pmax()<<"] Pa");
        }
        return true;
    }

    bool cSWEOSarg::CheckRanges_X(double *Xrange, double XMIN, double XMAX) {
        if(Xrange[0]>XMAX || Xrange[0]<XMIN)
        {
            CLI_FAIL("The minimum value of salinity specified by -R argument Xmin ="<<Xrange[0]<<", out of range ["<<XMIN<<", "<<XMAX<<"]");
        }
        if(Xrange[1]>XMAX || Xrange[1]<XMIN)
        {
            CLI_FAIL("The maximum value of salinity specified by -R argument Xmax ="<<Xrange[1]<<", out of range ["<<XMIN<<", "<<XMAX<<"]");
        }
        return true;
    }

                    bool cSWEOSarg::CheckEOS_Energy()
    {
        if (!(m_haveT || m_haveH))
        {
            CLI_USAGE("Selected calculation mode is "<<m_valueD<<"D calculation, change "<<m_valueV<<", but you didn't set a proper fixed temperature value by -T or set a proper fixed enthalpy by -H");
        }
        m_TorH = (m_haveT ? LOOKUPTABLE_FOREST::EOS_ENERGY_T : LOOKUPTABLE_FOREST::EOS_ENERGY_H);

        return true;
    }


}
