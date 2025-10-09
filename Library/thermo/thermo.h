/**
 * @file thermo.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Head file of xThermal
 * @version 0.1
 * @date 2022-01-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifndef XTHERMO_THERMO_H
#define XTHERMO_THERMO_H

#include "stdfunc.h"
#include <iostream>
#include <fstream>
#include "DataStructures.h"
#include "Exception.h"
#include <vector>
#include <cmath>
// it is fine for mac and linux with this head, but it is required in windows
#include <algorithm>
#include "triangle.h"
#include "LookUpTableForest.h"
#include <string>
#include "MultiProgressBar.h"

/**
 * @brief Namespace of xThermal library.
 * 
 */
namespace xThermal
{
    const double R = 461.51805;        /**< Specific gass constant: \f$ J/kg/K \f$ */
    // see also http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
    // map: fluid name to index
    // const std::map<std::string, int > FluidName2Index={"dsafds",2};
    // map: fluid index to fluid name
    enum INDEX_FLUID
    {
        Fluid_Unknown = -1,
        FLUID_Water,
        FLUID_H2O_NaCl
    };

    /**
     * @brief Top abstract class of the thermodynamic model.
     * 
     */
    class cxThermal
    {
    private:
        // MAP_PhaseRegion2Name m_map_phase2name;
        //LUT related data
        int m_num_threads;
        void *m_pLUT;
        int m_dim_lut;
        void *m_pLUT_lookup; //used as input LUT and used to lookup.
        int m_dim_lut_lookup;
        int m_index_Rho_in_LUT;
        std::map<int, propInfo> m_supported_props;
        std::map<int, propInfo> m_update_which_props;
        // LUT related methods
        void initialize_data();
        void init_supported_props();
        void parse_update_which_props(int update_which_props);
        void destroyLUT(void* pLUT, int& dim_lut);
        template<int dim>
        void interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf, double* xyz_min_target, ThermodynamicProperties& prop, const double xyz[dim]);

        template<int dim>
        void interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf, double* xyz_min_target, double* props, const double xyz[dim]);

    public:
        cxThermal(/* args */);
        virtual ~cxThermal();
        bool m_isShowProgressBar;
        void showProgressBar(bool isShow){m_isShowProgressBar = isShow;};
        static INDEX_FLUID validateFluid(const std::string& fluidName);
    public:
        virtual std::string name() = 0;             /**< Name of the model */
        // -------- Some required information for every child class
        virtual double Tmin() = 0;                  /**< Get the minimum temperature in K */
        virtual double Tmax() = 0;                  /**< Get the maximum temperature in K */
        virtual double pmin() = 0;                  /**< Get the minimum pressure in Pa */
        virtual double pmax() = 0;                  /**< Get the maximum pressure in Pa */
        virtual double Ttriple() = 0;               /**< Get the triple point temperature in K */
        virtual double T_critical() = 0;            /**< Return the critical temperature in K */
        virtual double p_critical() = 0;            /**< Return the critical pressure in Pa */
        virtual double rhomass_critical() = 0;      /**< Return the critical mass density in \f$ kg/m^3 \f$ */
        virtual double molar_mass() = 0;            /**< Return the molar mass in kg/mol */
        virtual double rhomolar_critical(){return rhomass_critical()/molar_mass();};     /**< Calculate the critical molar density in \f$ mol/m^3 \f$ */
        // Update 
        /**
         * @brief Update state using T, p as independent variables.
         * 
         * @param T [K]
         * @param p [Pa]
         */
        virtual void UpdateState_TPX(ThermodynamicProperties& props, const double& T, const double& p, const double& X=0) {throw NotImplementedError(name() + " does not implement UpdateState_TPX function");};
        virtual void UpdateState_HPX(ThermodynamicProperties& props, const double& H, const double& p, const double& X=0) {throw NotImplementedError(name() + " does not implement UpdateState_HPX function");};
        virtual ThermodynamicPropertiesVector UpdateState_TPX(const std::vector<double>& T, const std::vector<double>& p, const std::vector<double>& X, bool isMeshGrid=false) final;
        virtual ThermodynamicPropertiesVector UpdateState_HPX(const std::vector<double>& H, const std::vector<double>& p, const std::vector<double>& X, bool isMeshGrid=false) final;
        // for matlab API
        virtual void UpdateState_TPX(ThermodynamicPropertiesArray& stateArray, const size_t& N, const double* T, const double* p, const double* X=NULL);
        virtual void UpdateState_HPX(ThermodynamicPropertiesArray& stateArray, const size_t& N, const double* H, const double* p, const double* X=NULL);

        virtual std::string phase_name(PhaseRegion phase_index) final;

        // ====== LUT related functions ======
        virtual void set_num_threads(int num_threads) final;
        virtual int get_num_threads() final;
        virtual PhaseRegion findPhaseRegion_TPX(const double& T, const double& p, const double& X=0){throw NotImplementedError(name() + " does not implement findPhaseRegion_pTX function");};
        virtual const std::map<int, propInfo>& get_UpdateWhichProps() final;
        virtual void createLUT_2D(double xy_min[2], double xy_max[2], double constZ, LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var, LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
                                  int min_level = 4, int max_level = 6, int update_which_props=0) final;
        virtual void createLUT_2D(double xmin, double xmax, double ymin, double ymax, double constZ, LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var, LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
                                  int min_level = 4, int max_level = 6, int update_which_props=0)final;
        virtual void createLUT_3D(double xyz_min[3], double xyz_max[3], LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0) final;
        virtual void createLUT_3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0) final;
        virtual void save_lut_to_vtk(std::string filename, bool isNormalizeXYZ=true) final;
        virtual void save_lut_to_binary(std::string filename) final;
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *lookup(ThermodynamicProperties& prop, double x, double y);
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *lookup(double* props, double* xyz_min_target, double x, double y, bool is_cal=true);
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *lookup(ThermodynamicProperties& prop, double x, double y, double z);
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *lookup(double* props, double* xyz_min_target, double x, double y, double z, bool is_cal=true);
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *lookup_only(ThermodynamicProperties& prop, double x, double y);
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *lookup_only(ThermodynamicProperties& prop, double x, double y, double z);
        void loadLUT(const std::string& filename, bool printStatus=true);
        void* get_pLUT();
        void* get_pLUT_lookup();
        int get_dim_lut(){return m_dim_lut; };
        int get_dim_lut_lookup(){return m_dim_lut_lookup; };
        PhaseRegion Rho_lookup(double& Rho_estimate, double& Rho_min, double& Rho_max, const double& T, const double& P);
        //====================================
        /**
         * Calculate boiling pressure [Pa] of water for a given T [K]
         * @param T
         * @return
         */
        virtual double Boiling_p(const double& T) {throw NotImplementedError(name() + " does not implement Boiling_p(const double& T) function");};
        /**
         * Calculate boiling temperature [K] of water for a given p [Pa]
         * @param p
         * @return
         */
        virtual double Boiling_T(const double& p){throw NotImplementedError(name() + " does not implement Boiling_T(const double& p) function");};
        /**
         * Calculate both boiling temperature and density of liquid and vapor phase.
         * @param T [K]
         * @param rho_l [kg/m3]
         * @param rho_v [kg/m3]
         * @return
         */
        virtual double Boiling_p(const double& T, double & rho_l, double & rho_v){throw NotImplementedError(name() + " does not implement Boiling_p(const double& T, double & rho_l, double & rho_v) function");};
        /**
         * Calculate both boiling temperature and density of liquid and vapor phase.
         * @param p [Pa]
         * @param rho_l [kg/m3]
         * @param rho_v [kg/m3]
         * @return
         */
        virtual double Boiling_T(const double& p, double& rho_l, double& rho_v){throw NotImplementedError(name() + " does not implement Boiling_T(const double& p, double& rho_l, double& rho_v) function");};
        virtual double Boiling_T(const double& p, ThermodynamicProperties& props){throw NotImplementedError(name() + " does not implement Boiling_T(const double& p, ThermodynamicProperties& props) function");};
        virtual double Boiling_p(const double& T, ThermodynamicProperties& props){throw NotImplementedError(name() + " does not implement Boiling_p(const double& T, ThermodynamicProperties& props) function");};
        void Triangulation(const std::vector<double>& x_poly, const std::vector<double>& y_poly, const double pointInMesh[2], const double dxdy[2], TriMesh& trimesh, bool normalize=true);
        virtual void writeTriMesh2Txt(const TriMesh& mesh, std::string path_out=".")final;
        // write mesh to vtu file
        virtual void writeXXYYZZ2VTU(std::string vtuFile, const std::vector<std::vector<double> >& XX, const std::vector<std::vector<double> >& YY, const std::vector<std::vector<double> >& ZZ,
                const double scale_x=1.0, const double scale_y=1.0, double scale_z=1.0)final;
        virtual void writeLine2VTU(std::string vtuFile, const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z, const double scale_x=1.0, const double scale_y=1.0, double scale_z=1.0)final;
        virtual void writePhaseBoundaries2VTU(std::string outputPath, const PhaseBoundaries& phaseBoundaries, double scale_T=1, double scale_p = 1, double scale_X = 1)final;
        virtual void normalizePhaseBoundaries(PhaseBoundaries& phaseBoundaries)final;
        virtual void writeMeshGrid2VTK(const std::string &vtkFile, const std::vector<double> &x, const std::string& xTitle, const std::vector<double> &y, const std::string& yTitle, const std::vector<double> &z, const std::string& zTitle,
                                       const std::vector<std::vector<double>> &props, const std::vector<propInfo> &propsInfo, bool isNormalize=false) final;
    public:
        template<class T> T min_vector(const std::vector<T> &x)
        {
            T min = 1e40;
            std::size_t N = x.size();
            for (std::size_t i = 0; i < N; ++i)
            {
                T axi = x[i];
                if (axi < min){ min = axi; }
            }
            return min;
        };
        template<class T> T max_vector(const std::vector<T> &x)
        {
            T max = -1E30;
            std::size_t N = x.size();
            for (std::size_t i = 0; i < N; ++i)
            {
                T axi = x[i];
                if (axi > max){ max = axi; }
            }
            return max;
        };
        template<class T> T sum_vector(const std::vector<T> &x)
        {
            T sum = 0;
            std::size_t N = x.size();
            for (std::size_t i = 0; i < N; ++i)
            {
                sum += x[i];
            }
            return sum;
        };
        template<class T> T mean_vector(const std::vector<T> &x)
        {
            return sum_vector(x)/x.size();
        };
        /// Make a linearly spaced vector of points
        template <typename T> std::vector<T> linspace(T xmin, T xmax, std::size_t n, bool isLogScale=false) {
            std::vector<T> x(n, 0.0);
            // if one of xmin and xmax is negative, logscale is not available
            if(xmin<0 || xmax<0)isLogScale=false;
            if(isLogScale)
            {
                double delta = (log10(xmax)-log10(xmin))/(n-1);
                double min = log10(xmin);
                for ( std::size_t i = 0;  i < n; ++i) {
                    x[i] = pow(10.0, delta*i+min);
                }
            }else
            {
                for ( std::size_t i = 0;  i < n; ++i) {
                    x[i] = (xmax-xmin)/(n-1)*i+xmin;
                }
            }
            return x;
        }
        template <typename T> std::vector<T> linspace(T xmindxxmax[3]) {
            double dx = xmindxxmax[1];
            size_t N = (size_t)((xmindxxmax[2] - xmindxxmax[0])/dx + 1);
            double dx_new = (xmindxxmax[2] - xmindxxmax[0])/(N -1);
            std::vector<T> x(N);
            for (size_t i = 0; i < N; ++i) {
                x[i] = xmindxxmax[0] + dx_new*i;
            }
            // for (T x0 = xmindxxmax[0]; x0 <=xmindxxmax[2]; x0 +=xmindxxmax[1]) {
            //     x.push_back(x0);
            // }
            return x;
        }
        /**
        * Create a regular mesh grid.
        * @tparam T Data type
        * @param xmindxxmax [xmin, dx, xmax]
        * @param ymindyymax [ymin, dy, ymax]
        * @param XX 1D vector of X in order of [x0, x1, ..., xn; x0, x1, ..., xn; ... ; x0, x1, ..., xn]
        * @param YY 1D vector of Y in order of [y0, y0, ..., y0; y1, y1, ..., y1; ...; yn, yn, ..., yn]
        * @return [nx, ny]
        */
        template <typename T> std::vector<size_t> meshgrid(T xmindxxmax[3], T ymindyymax[3], std::vector<T>& XX, std::vector<T>& YY) {
            std::vector<T> x, y;
            for (T x0 = xmindxxmax[0]; x0 <=xmindxxmax[2]; x0 +=xmindxxmax[1]) x.push_back(x0);
            for (T y0 = ymindyymax[0]; y0 <=ymindyymax[2]; y0 +=ymindyymax[1]) y.push_back(y0);
            size_t  nx = x.size(), ny = y.size();
            XX.resize(nx * ny);
            YY.resize(XX.size());
            // make grid
            int ind = 0;
            for (size_t i = 0; i < ny; ++i) {
                for (size_t j = 0; j < nx; ++j) {
                    ind = j + i*nx;
                    XX[ind] = x[j];
                    YY[ind] = y[i];
                }
            }
            return {nx, ny};
        }

    public:  //for python wrapper
        ThermodynamicProperties Boiling_p_props(const double& T);
        ThermodynamicProperties Boiling_T_props(const double& p);
        ThermodynamicProperties UpdateState_TPX(const double& T, const double& p, const double& X=0);
        ThermodynamicProperties UpdateState_HPX(const double& H, const double& p, const double& X=0);
        Head_AMR_LUT getLutInfo(std::string file_lut, bool printSummary=true);
    };

    //LUT related function
    void fill_prop2data(cxThermal* pEOS, const ThermodynamicProperties* prop, const std::map<int, propInfo>& update_which_props, double* data);
}

#endif
