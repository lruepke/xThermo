

# File LookUpTableForest.h

[**File List**](files.md) **>** [**AMR\_LUT**](dir_eb711624b693511de215f59cb195e6f2.md) **>** [**LookUpTableForest.h**](LookUpTableForest_8h.md)

[Go to the documentation of this file](LookUpTableForest_8h.md)


```C++
#ifndef LOOKUPTABLEFOREST_H
#define LOOKUPTABLEFOREST_H
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "stdfunc.h"
#include <cmath>

#if USE_OMP == 1
    #include <omp.h>
#endif

#include "thermo.h"

namespace LOOKUPTABLE_FOREST
{
    #define MAX_FOREST_LEVEL 29
    #define ExtensionName_PointIndexFile "pi"
    struct Quad_index
    {
        int i = 0, j = 0, k = 0;
        bool operator< (const Quad_index &ijk) const
        {
            return i < ijk.i || (i==ijk.i && j<ijk.j) || (i==ijk.i && j==ijk.j && k<ijk.k);
        }
    };
    typedef unsigned int int_pointIndex;
    struct PropsData
    {
        double** data = NULL;
        int_pointIndex num_points = 0;
        int num_props = 0;
        void create()
        {
            if(num_props>0)
            {
                data = new double*[num_points];
                for (int_pointIndex i = 0; i < num_points; i++)
                {
                    data[i] = new double[num_props];
                }
            }
        }
        void clear()
        {
            if(data) //need to check where the data array is created or not, e.g., lutInfo, doesn't create data array.
            {
                for (int_pointIndex i = 0; i < num_points; i++)
                {
                    delete[] data[i];
                    data[i] = NULL;
                }
                delete[] data;
                num_points = 0;
            }
        }
    };
    
    
    // non complete declaration
    template <int dim, typename USER_DATA> struct NonLeafQuad;
    template <int dim, typename USER_DATA> struct LeafQuad;
    template <int dim, typename USER_DATA>union Quad_data
    {
        NonLeafQuad<dim, USER_DATA>* nonleaf;
        LeafQuad<dim, USER_DATA>* leaf = NULL;
    };
    template <int dim, typename USER_DATA>
    struct Quadrant
    {
        unsigned char       level;
        bool                isHasChildren;
        Quad_data<dim, USER_DATA> qData; //leaf quad and nonleaf quad use different type of data
    };

    // data for non-leaf quad
    template <int dim, typename USER_DATA> 
    struct NonLeafQuad
    {
        Quadrant<dim, USER_DATA>   *children[1<<dim];
    };
    // data for leaf quad
    template <int dim, typename USER_DATA> 
    struct LeafQuad
    {
        Quadrant<dim, USER_DATA>* parent =NULL;
        USER_DATA           *user_data = NULL;
        unsigned int index_props[1<<dim]; //index of property on each node.
    };
    
    enum CONST_WHICH_VAR 
    {
        CONST_NO_VAR_TorHPX,     
        CONST_TorH_VAR_XP,      
        CONST_P_VAR_XTorH,      
        CONST_X_VAR_TorHP       
    }; //only used for 2D case, CONST_NO means 3D
    
    enum NeedRefine {NeedRefine_NoNeed, NeedRefine_PhaseBoundary, NeedRefine_Rho, NeedRefine_H, NeedRefine_Mu};
    struct RMSD_RefineCriterion
    {
        double Rho;
        double H;
        double Mu;
    };

    // add data struct definition for different system, e.g., H2ONaCl. Actually we can move this data type definition to H2ONaCl.H, because LookUpTableForest class never care about this data type definite, it just accept whatever data type through template argument. But for the TCL API, if move this to other place, it will cause some compiling errors. So keep it here before finding better solution.
    
    enum EOS_ENERGY {
        EOS_ENERGY_T, 
        EOS_ENERGY_H  
        };
    template <int dim>
    struct FIELD_DATA
    {
        // point data field
        // H2ONaCl::PROP_H2ONaCl prop_point[1<<dim]; // properties at vertiex
        // cell data field
        NeedRefine need_refine; // indicator of what kind of the need-refined quad position (phase boundary), or what kind of properties need to refine
        xThermal::PhaseRegion phaseRegion_cell;
    };

    inline int get_dim_from_binary(std::string filename)
    {
        FILE* fpin = NULL;
        fpin = fopen(filename.c_str(), "rb");
        if(!fpin)ERROR("Open file failed: "+filename);
        int dim0;
        fread(&dim0, sizeof(dim0), 1, fpin);
        fclose(fpin);
        return dim0;
    }
    template <int dim, typename USER_DATA> 
    class LookUpTableForest
    {
    private:
        long int m_num_quads;
        int m_num_leaves;
        int m_num_need_refine;
        size_t  m_data_size;
        double m_length_scale[dim]; 
        Quadrant<dim,USER_DATA> m_root;
        void init_Root(Quadrant<dim,USER_DATA>& quad);
        void release_quadrant_data(Quadrant<dim,USER_DATA>* quad);
        void release_leaves(Quadrant<dim,USER_DATA>* quad);
        void getLeaves(std::vector<Quadrant<dim,USER_DATA>* >& leaves, long int& quad_counts, Quadrant<dim,USER_DATA>* quad);
        void refine(Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level));
        void write_vtk_cellData(std::ofstream* fout, std::string type, std::string name, std::string format);
        void searchQuadrant(Quadrant<dim,USER_DATA>* quad_source, Quadrant<dim,USER_DATA> *&quad_target, double* xyz_min_target, double x_ref, double y_ref, double z_ref);
        void init(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer);
        void write_forest(FILE* fpout, FILE* fpout_point_index, Quadrant<dim,USER_DATA>* quad, int order_child, bool is_write_data);
        void read_forest(FILE* fpin_forest, FILE* fpin_point_index, Quadrant<dim,USER_DATA>* quad, int order_child);
        void construct_map2dat();
        void get_unique_points_leaves(std::map<Quad_index, int_pointIndex>& map_unique_points, int& num_leaves, long int& num_quads, int& num_need_refine, Quadrant<dim,USER_DATA>* quad, Quad_index ijk_quad, unsigned int length_quad);
        void pass_props_pointer_leaves(std::map<Quad_index, int_pointIndex>& map_unique_points, Quadrant<dim,USER_DATA>* quad, Quad_index ijk_quad, unsigned int length_quad);
        void read_props_from_binary(const std::string& filename_forest, bool printStatus=true);
        bool read_forest_from_binary(const std::string& filename, bool read_only_header=false, bool printStatus=true);
    public:
        void    *m_eosPointer;      //pass pointer of EOS object (e.g., the pointer of a object of cH2ONaCl class) to the forest through construct function, this will give access of EOS stuff in the refine call back function, e.g., calculate phase index and properties
        double  m_constZ;         // only valid when dim==2, i.e., 2D case, the constant value of third dimension, e.g. in T-P space with constant X.
        int     m_min_level;
        int     m_max_level;
        double m_xyz_min[dim];
        double m_xyz_max[dim];
        int     m_num_children;
        int     m_num_node_per_quad; //how many nodes will be stored in a quad: only for data storage
        std::map<int, xThermal::propInfo> m_map_props;
        std::map<int, int> m_map_prop2index; // map property name macro(e.g., Update_prop_rho) to index, i.e. index of this property in m_props_unique_points_leaves.data[indNode][index of prop]
        PropsData m_props_unique_points_leaves;
        CONST_WHICH_VAR m_const_which_var; 
        EOS_ENERGY m_TorH; 
        // double  m_physical_length_quad[MAX_FOREST_LEVEL][dim]; //Optimization: store the length of quad in each dimension as a member data of the forest, therefore don't need to calculate length of quad, just access this 2D array according to the quad level. 
        RMSD_RefineCriterion m_RMSD_RefineCriterion;
        inline void set_min_level(int min_level){m_min_level = min_level;};
        Quadrant<dim,USER_DATA>* get_root(){return &m_root;};
        // int searchQuadrant(double x, double y, double z);
        void searchQuadrant(Quadrant<dim,USER_DATA> *&targetLeaf, double* xyz_min_target, double x, double y, double z);
        void get_quadrant_physical_length(int level, double physical_length[dim]);
        void refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level));
        void get_ijk_nodes_quadrant(Quadrant<dim,USER_DATA>* quad, const Quad_index* ijk_quad, int num_nodes_per_quad, Quad_index* ijk);
        void assemble_data(void (*cal_prop)(LookUpTableForest<dim,USER_DATA>* forest, std::map<Quad_index, double*>& map_ijk2data));
        void construct_props_leaves(void (*cal_prop)(LookUpTableForest<dim,USER_DATA>* forest, std::map<Quad_index, unsigned int>& map_ijk2data, double** data));
        void ijk2xyz(const Quad_index* ijk, double& x, double& y, double& z);
        void union_ijk2xyz(Quadrant<dim, USER_DATA>* quad, Quad_index& ijk_backup);
        void write_to_vtk(std::string filename, bool write_data=true, bool isNormalizeXYZ=true);
        void write_to_binary(std::string filename, bool is_write_data=true);
        void write_point_index(std::string filename_forest);
        void write_point_index(FILE* fpout_point_index, Quadrant<dim,USER_DATA>* quad);
        void print_summary();
        int get_num_leaves(){return m_num_leaves;};
        int get_num_quads(){return m_num_quads;};
        int get_num_need_refine(){return m_num_need_refine;};
        std::string byte2string(double bytes);
        LookUpTableForest(double xyz_min[dim], double xyz_max[dim], EOS_ENERGY TorH, int max_level, std::map<int, xThermal::propInfo> name_props, void* eosPointer=NULL); //3D case
        LookUpTableForest(double xy_min[dim], double xy_max[dim], double constZ, CONST_WHICH_VAR const_which_var, EOS_ENERGY TorH, int max_level, std::map<int, xThermal::propInfo> name_props, void* eosPointer=NULL); //2D case
        LookUpTableForest(std::string filename_forest, void* pointer=NULL, bool printStatus=true); //load from exist binary file
        void destory();
        ~LookUpTableForest();
    };

    typedef LookUpTableForest<2, FIELD_DATA<2> > LookUpTableForest_2D;
    typedef LookUpTableForest<3, FIELD_DATA<3> > LookUpTableForest_3D;
    #define CHECK_REFINE_PROP_RMSD(PROP) \
    { \
        double mean_Prop = props_refine_check[0].PROP; \
        for(int i=0;i<forest->m_num_children;i++)mean_Prop += props_refine_check[i+1].PROP; \
        mean_Prop = mean_Prop / (forest->m_num_children + 1);  \
        double RMSD_Prop = pow((props_refine_check[0].PROP - mean_Prop)/mean_Prop, 2.0); \
        for(int i=0;i<forest->m_num_children;i++)RMSD_Prop += pow((props_refine_check[i+1].PROP - mean_Prop)/mean_Prop, 2.0); \
        RMSD_Prop = sqrt(RMSD_Prop/(forest->m_num_children + 1)); \
        if(RMSD_Prop > forest->m_RMSD_RefineCriterion.PROP) \
        { \
            data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_##PROP; \
        } \
    }
}


#endif
```


