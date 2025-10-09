%module H2ONaCl
%{
    // #define SWIG_FILE_WITH_INIT

    #include "stdfunc.h"
    #include "thermo.h"
    #include "H2ONaCl.h"
    #include <vector>
    #include <algorithm>
    #include <sstream>
    using namespace xThermal;
%}

%include std_string.i
%include std_vector.i
%include std_map.i
%include <std_string.i>
// convert std::Vector to numpy array
// Cmake option: USE_NUMPY
#ifdef HAVE_numpy
//https://github.com/numpy/numpy/blob/main/tools/swig/numpy.i
%include "numpy.i"
%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};
#endif

%inline %{
using namespace std;
%}
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<string>;
%template(ConstCharVector) std::vector<const char*>;
%template(UnsignedLongVector) std::vector<unsigned long int>;
%template(vector_vector_double) std::vector<std::vector<double> >;
%template(vector_vector_int) std::vector<std::vector<int> >;
%template(map_string_string) std::map<std::string, std::string>;
// type map
%include "typemaps.i"
//see https://github.com/swig/swig/blob/master/Examples/test-suite/li_std_vector.i for more templates

enum PhaseRegion
{
    MixPhaseRegion=-1,      /**< MixPhaseRegion is used in AMR lookup table. */
    SinglePhase_L,          /**< Single phase: liquid. Regard single-phase "liquid-like" and "vapor-like" fluid as "Single phase liquid" as well. */
    SinglePhase_V,          /**< Single phase: vapor*/
    SinglePhase_S,          /**< Single phase: solid, e.g. NaCl */
    Supercritical,          /**< Supercritical (p > pc, T > Tc) for pure water */
    Supercritical_vapor,    /**< Supercritical vapor (p < pc, T > Tc) for pure water */
    Supercritical_liquid,   /**< Supercritical liquid (p > pc, T < Tc) for pure water */
    Critical,               /**< At the critical point of pure water or on critical curve of H2ONaCl */
    TwoPhase_VL_Water,           /**< Two phase: liquid + vapor, salinity = 0*/
    TwoPhase_LH,           /**< Two phase: liquid + halite */
    TwoPhase_VH,           /**< Two phase: vapor + halite */
    // TwoPhase_V_L_L,         /**< Two phase: vapor + liquid coexistence, liquid branch */
    // TwoPhase_V_L_V,         /**< Two phase: vapor + liquid coexistence, vapor branch */
    TwoPhase_VL,            /**< Two phase: VL region */
    ThreePhase_VLH,       /**< Three phase: vapor + liquid + halite */
    Unknown,                /**< Unknown phase region for bug report*/
    NotImposed              /**< Phase is not imposed */
};
struct ThermodynamicProperties
{
    double T=0, p=0, X=0, H=0;
    PhaseRegion phase=Unknown;
    double S_l=0, S_v=0, S_h=0;
    double X_l=0, X_v=0;
    double Rho_l=0, Rho_v=0, Rho_h=0;
    double H_l=0, H_v=0, H_h=0;
    double Cp_l = 0, Cp_v = 0, Cp_h = 0;
    double Mu_l = 1, Mu_v = 1;
    double Rho=0, Cp = 0, Mu = 1;
    // derivative
    double dRhodP = 0, dRhodP_l = 0, dRhodP_v = 0;
    double dRhodT = 0, dRhodT_l = 0, dRhodT_v = 0;
    // commonly used derivatives for hydrothermal convection simulation
    double IsothermalCompressibility = 0, IsothermalCompressibility_l =0, IsothermalCompressibility_v =0;
    double IsobaricExpansivity=0, IsobaricExpansivity_l =0, IsobaricExpansivity_v =0;
    std::string fluidName;
    std::string info(const ThermodynamicProperties& self);
};
%extend ThermodynamicProperties {
        std::string __str__() {
            return self->info(*self);
        }
};
struct ThermodynamicPropertiesVector
{
    std::vector<double> T, p, X, H;
    std::vector<int> phase;
    std::vector<double> S_l, S_v, S_h;
    std::vector<double> X_l, X_v;
    std::vector<double> Rho_l, Rho_v, Rho_h;
    std::vector<double> H_l, H_v, H_h;
    std::vector<double> Cp_l, Cp_v, Cp_h;
    std::vector<double> Mu_l, Mu_v;
    std::vector<double> Rho, Cp, Mu;
    std::vector<double> dRhodP, dRhodP_l, dRhodP_v;
    std::vector<double> dRhodT, dRhodT_l, dRhodT_v;
    std::vector<double> IsothermalCompressibility, IsothermalCompressibility_l, IsothermalCompressibility_v;
    std::vector<double> IsobaricExpansivity, IsobaricExpansivity_l, IsobaricExpansivity_v;
    std::string fluidName;
    std::string info(const ThermodynamicPropertiesVector& self);
};
%extend ThermodynamicPropertiesVector {
        std::string __str__() {
            return self->info(*self);
        }
};

enum PhaseType{Vapor, Liquid};
/**
     * 2D Deformable linear "structured" mesh. Related to matplotlib ax.plot_wireframe()
     */
struct DeformLinearMesh
{
    std::vector<std::vector<double> > T, p, X;
};

/**
 * 2D triangle mesh structure. Relate to matplotlib ax.tri_plot() in 2D or ax.plot_trisurf() in 3D
 */
struct TriMesh
{
    std::vector<double> x,y,z;
    std::vector<std::vector<int> > connection;
    std::vector<std::string>names={"x","y","z"};
};
/**
* Data structure for a phase boundary surface
*/
struct Surface
{
    std::string name, shortName;
    DeformLinearMesh mesh;
    std::string color;
};
/**
* Data structure for a phase boundary line
*/
struct Line
{
    std::string name,shortName;
    std::vector<double > T,p,X;
    std::string color = "black";
};
/**
* Data structure for a phase boundary point
*/
struct Point
{
    std::string name, shortName;
    double T,p,X;
    std::string color = "black";
};
#define SCALE_X_linear 1
#define SCALE_X_log 2
#define SCALE_X_loglinear 3
struct PhaseBoundaries
{
    std::vector<Surface> surfaces;
    std::vector<Line> lines;
    std::vector<Point> points;
    int scale_X; //"linear", "log", "loglinear"
    double Xcenter; //The log and linear part is seperated by Xcenter: kg/kg
    double ratio_log_to_linear;
    double Tmin, Tmax, pmin,pmax, Xmin, Xmax;
};

/**
 * Head info of AMR_LUT which are used for Python API
 */
struct Head_AMR_LUT
{
    int dim=0;
    int space=-1;
    std::string spaceName="Unknown";
    int constWhich=-1;
    std::string constWhich_Name="Unknown";
    double constValue=0; //SI unit
    std::vector<std::string> var_names;
    std::vector<std::vector<double> > var_ranges;//SI unit
    std::vector<double> var_maxResolution; //SI unit
    int min_level, max_level;
    int num_leaves, num_quads, num_nodes, num_need_refine;
    int num_leaves_nextRefine;
    int num_props;
    std::vector<std::string> prop_names;
    std::string memory_total, memory_leaves, memory_nonLeaves, memory_quads, memory_properties;
};

struct COLOR
{
    double r,g,b;
};
// 这个naturalvar的标注非常重要，否则在python中无法访问结构体里面的vector，比如T,P,X, H
%naturalvar Polygon_slice::T;
%naturalvar Polygon_slice::P;
%naturalvar Polygon_slice::X;
%naturalvar Polygon_slice::H;
%naturalvar Polygon_slice::fc;
%naturalvar Polygon_slice::ec;
struct Polygon_slice
{
    std::vector<double> T,P,X,H;
    std::string name;
    std::vector<double> fc; //face color: [R,G,B]
    std::vector<double> ec; //edge color: [R,G,B]
    Polygon_slice(){}
    Polygon_slice(std::vector<double> T0, std::vector<double> P0, std::vector<double> X0, std::vector<double> H0, std::string name0, COLOR fc0, COLOR ec0)
    : T(std::move(T0)), P(std::move(P0)), X(std::move(X0)), H(std::move(H0)), name(std::move(name0)), fc{fc0.r, fc0.g, fc0.b}, ec{ec0.r, ec0.g, ec0.b}
    {

    }
    bool operator<(const Polygon_slice &other) const { return T < other.T; }
    bool operator==(const Polygon_slice &other) const { return name == other.name; }
    std::string info(const Polygon_slice& self);
};
%naturalvar Line_slice::T;
%naturalvar Line_slice::P;
%naturalvar Line_slice::X;
%naturalvar Line_slice::H;
%naturalvar Line_slice::color;
struct Line_slice
{
    std::vector<double>T,P,X,H;
    std::string name;
    std::vector<double> color;
    std::string ls;
    double lw;
    Line_slice(){}
    Line_slice(std::vector<double> T0, std::vector<double> P0, std::vector<double> X0, std::vector<double> H0, std::string name, COLOR c, std::string ls="solid", double lw=1)
    : T(std::move(T0)),P(std::move(P0)),X(std::move(X0)),H(std::move(H0)),name(std::move(name)),color{c.r, c.g, c.b}, ls(std::move(ls)),lw(lw)
    {
    }
    bool operator<(const Line_slice &other) const { return T < other.T; }
    bool operator==(const Line_slice &other) const { return name == other.name; }
    std::string info(const Line_slice& self);
};
%naturalvar Point_slice::mfc;
%naturalvar Point_slice::mec;
struct Point_slice
{
    double T,P,X,H;
    std::string name;
    std::string marker;
    std::vector<double> mfc;
    std::vector<double> mec;
    // the construct function is needed for Python API
    Point_slice(){}
    Point_slice(double T, double P, double X, double H, std::string name, std::string marker, COLOR mfc0, COLOR mec0)
            : T(T),P(P), X(X), H(H), name(std::move(name)), marker(std::move(marker)), mfc{mfc0.r, mfc0.g, mfc0.b}, mec{mec0.r, mec0.g, mec0.b}
    {
    }
    bool operator<(const Point_slice &other) const { return T < other.T; }
    bool operator==(const Point_slice &other) const { return name == other.name; }
    std::string info(const Point_slice& self);
};
%extend Point_slice {
        std::string __str__() {
            return self->info(*self);
        }
};
%extend Line_slice {
        std::string __str__() {
            return self->info(*self);
        }
};
%extend Polygon_slice {
        std::string __str__() {
            return self->info(*self);
        }
};
%template(Point_sliceVector) std::vector<Point_slice>;
%template(Line_sliceVector) std::vector<Line_slice>;
%template(Polygon_sliceVector) std::vector<Polygon_slice>;
%template(stringPolygon_sliceMap) std::map<string, Polygon_slice>;
%template(stringLine_sliceMap) std::map<string, Line_slice>;
%template(stringPoint_sliceMap) std::map<string, Point_slice>;
%template(stringPoint_sliceVectorMap) std::map<string, std::vector<Point_slice>>;
%template(stringLine_sliceVectorMap) std::map<string, std::vector<Line_slice>>;
%template(stringPolygon_sliceVectorMap) std::map<string, std::vector<Polygon_slice>>;

//extend of vector
%extend std::vector<Point_slice> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"Vector of <class>Point_slice, size="<<self->size()<<std::endl;
            return sin.str();
        }
};
%extend std::vector<Line_slice> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"Vector of <class>Line_slice, size="<<self->size()<<std::endl;
            return sin.str();
        }
};
%extend std::vector<Polygon_slice> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"Vector of <class>Polygon_slice, size="<<self->size()<<std::endl;
            return sin.str();
        }
};
//extend of map vector
%extend std::map<string, std::vector<Point_slice>> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"{";
            for(auto &item : *self)sin<<"'"<<item.first<<"' ";
            sin<<"}"<<std::endl;
            return sin.str();
        }
};
%extend std::map<string, std::vector<Line_slice>> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"{";
            for(auto &item : *self)sin<<"'"<<item.first<<"' ";
            sin<<"}"<<std::endl;
            return sin.str();
        }
};
%extend std::map<string, std::vector<Polygon_slice>> {
        std::string __str__() {
            std::stringstream sin;
            sin<<"{";
            for(auto &item : *self)sin<<"'"<<item.first<<"' ";
            sin<<"}"<<std::endl;
            return sin.str();
        }
};
/**
 * Struct of phase region of a slice at constant T, P, H, or X. Contains phase region polygon, and some visualization properties.
 */
struct PhaseRegion_Slice
{
    std::map<std::string, std::vector<Polygon_slice> > regions;
    std::map<std::string, std::vector<Line_slice> > lines;
    std::map<std::string, std::vector<Point_slice> > points;
    std::string info(const PhaseRegion_Slice& slice);
};
%extend PhaseRegion_Slice {
    std::string __str__() {
            return self->info(*self);
        }
};

// ========================================================================
%apply double *OUTPUT {double& T_crit};
%apply double *OUTPUT {double& P_crit};
%apply double *OUTPUT {double& X_crit};
%apply double *OUTPUT {double& T_min, double& T_max};
%apply std::vector<double> *OUTPUT {std::vector<double>& P_crit, std::vector<double>& X_crit};
%apply std::vector<double> *OUTPUT {std::vector<double>& T_crit, std::vector<double>& X_crit};
%apply double *OUTPUT {double& P_crit, double& X_crit};
%apply double *OUTPUT {double& T_crit, double& X_crit};
%apply std::vector<double> *OUTPUT {std::vector<double>& T_crit};
%apply std::vector<double> *OUTPUT {std::vector<double>& P_crit};
%apply std::vector<double> *OUTPUT {std::vector<double>& res};
%apply double *OUTPUT {double& Tmin, double& Tmax};
%apply double *OUTPUT {double& X_l, double& X_v};
%apply std::vector<double> *OUTPUT {std::vector<double>& X_l, std::vector<double>& X_v};
%apply double *OUTPUT {double& X_V, double& X_L};
%apply std::vector<double> *OUTPUT {std::vector<double>& X_V, std::vector<double>& X_L};
%apply std::vector<double> *OUTPUT {std::vector<double>& n1, std::vector<double>& n2};
%apply std::vector<double> *OUTPUT {std::vector<double>& q1, std::vector<double>& q2};
%apply std::vector<double> *OUTPUT {std::vector<double>& P};
%apply double *OUTPUT {double& Hmin0, double& Hmax0};

//ignore some functions
//%ignore Boiling_p(const double& T, double&, double&);
//%ignore Boiling_T(const double& p, double&, double&);
%ignore UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);
%ignore UpdateState_HPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0);
////rename overloaded function
//%rename(Boiling_p_rho) Boiling_p(const double& T, double&, double&);
//%rename(Boiling_T_rho) Boiling_T(const double& p, double&, double&);


//ignore some warnings
%warnfilter(560,562) ;

#define xTHERMO_VAR
%include "thermo.h"
%include "H2ONaCl.h"