/**
 * @file H2ONaCl_LUT_RefineFuncI.H
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief In order to easily manage the code file, collect all the refinement function for the H2ONaCl lookup table(LUT) to this file.
 * These refinement function can not be implement in the .cpp file, so just include this file in \file H2ONaCl.cpp
 * @version 0.1
 * @date 2021-12-03
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef H2ONACL_LUT_REFINEFUNCI_H
#define H2ONACL_LUT_REFINEFUNCI_H

namespace xThermal
{

    /**
     * @brief Table refine function for the 2D case in T-P space with constant salinity.
     *
     * @tparam dim
     * @tparam USER_DATA
     * @param forest
     * @param quad
     * @param max_level
     * @return true
     * @return false
     */
    template <int dim, typename USER_DATA>
    bool RefineFunc_PTX(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, LOOKUPTABLE_FOREST::Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level)
    {
        // cout<<"level: "<<quad->level<<endl;
        if(quad->isHasChildren) return true; //if a quad has children, of course it need refine, but we don't need do anything at here, just return true.

        cxThermal* eosPointer =(cxThermal*)(forest->m_eosPointer); //read only! please DO NOT use this pointer to change any data in the EOS object!!! although it can change the member data.
        // double xv, xl;
        USER_DATA       *data = (USER_DATA *) quad->qData.leaf->user_data;
        data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_NoNeed; //initialize need_refine to be NoNeed
        // const double* physical_length_quad = forest->m_physical_length_quad[quad->level];
        double physical_length_quad[3];
        forest->get_quadrant_physical_length(quad->level, physical_length_quad);
        const int num_sample_x =2; //if want to make more safe check, could increase this number to generate more sample points, for most of cases, the default value 2 is enough.
        const int num_sample_y =2;
        const int num_sample_z = dim==3 ? 2 : 1;
        double dx_physical_quad = physical_length_quad[0] / (num_sample_x - 1.0);
        double dy_physical_quad = physical_length_quad[1] / (num_sample_y - 1.0);
        double x_physical_quad, y_physical_quad;
        PhaseRegion regionIndex[num_sample_x*num_sample_y*num_sample_z];
        ThermodynamicProperties* props_refine_check = new ThermodynamicProperties[forest->m_num_children + 1];
        switch (forest->m_const_which_var)
        {
            case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        regionIndex[iy*num_sample_x + ix] = eosPointer->findPhaseRegion_TPX(x_physical_quad, y_physical_quad, forest->m_constZ);
                    }
                }
                data->phaseRegion_cell     = eosPointer->findPhaseRegion_TPX(xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0, forest->m_constZ);

                // ========================================================
                // calculate properties: four vertices and one midpoint
                eosPointer->UpdateState_TPX(props_refine_check[0], xmin_quad,                               ymin_quad,                                forest->m_constZ);      //xmin,ymin
                eosPointer->UpdateState_TPX(props_refine_check[1], xmin_quad + physical_length_quad[0],     ymin_quad,                                forest->m_constZ);      //xmax,ymin
                eosPointer->UpdateState_TPX(props_refine_check[2], xmin_quad,                               ymin_quad + physical_length_quad[1],      forest->m_constZ);      //xmin,ymax
                eosPointer->UpdateState_TPX(props_refine_check[3], xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      forest->m_constZ);      //xmax,ymax
                eosPointer->UpdateState_TPX(props_refine_check[4], xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0,  forest->m_constZ);      //xc,yc
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        regionIndex[iy*num_sample_x + ix] = eosPointer->findPhaseRegion_TPX(forest->m_constZ, y_physical_quad, x_physical_quad);
                    }
                }
                data->phaseRegion_cell     = eosPointer->findPhaseRegion_TPX(forest->m_constZ, ymin_quad + physical_length_quad[1]/2.0, xmin_quad + physical_length_quad[0]/2.0);
                // ========================================================
                // calculate properties: four vertices and one midpoint
                eosPointer->UpdateState_TPX(props_refine_check[0], forest->m_constZ,   ymin_quad,                                xmin_quad);                                  //xmin,ymin
                eosPointer->UpdateState_TPX(props_refine_check[1], forest->m_constZ,   ymin_quad,                                xmin_quad + physical_length_quad[0]);        //xmax,ymin
                eosPointer->UpdateState_TPX(props_refine_check[2], forest->m_constZ,   ymin_quad + physical_length_quad[1],      xmin_quad);                                  //xmin,ymax
                eosPointer->UpdateState_TPX(props_refine_check[3], forest->m_constZ,   ymin_quad + physical_length_quad[1],      xmin_quad + physical_length_quad[0]);        //xmax,ymax
                eosPointer->UpdateState_TPX(props_refine_check[4], forest->m_constZ,   ymin_quad + physical_length_quad[1]/2.0,  xmin_quad + physical_length_quad[0]/2.0);    //xc,yc

            }
                break;
            case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        regionIndex[iy*num_sample_x + ix] = eosPointer->findPhaseRegion_TPX(y_physical_quad, forest->m_constZ, x_physical_quad);
                    }
                }
                data->phaseRegion_cell     = eosPointer->findPhaseRegion_TPX(ymin_quad + physical_length_quad[1]/2.0, forest->m_constZ, xmin_quad + physical_length_quad[0]/2.0);
                // ========================================================
                // calculate properties: four vertices and one midpoint
                eosPointer->UpdateState_TPX(props_refine_check[0], ymin_quad,                               forest->m_constZ, xmin_quad);                                  //xmin,ymin
                eosPointer->UpdateState_TPX(props_refine_check[1], ymin_quad,                               forest->m_constZ, xmin_quad + physical_length_quad[0]);        //xmax,ymin
                eosPointer->UpdateState_TPX(props_refine_check[2], ymin_quad + physical_length_quad[1],     forest->m_constZ, xmin_quad);                                  //xmin,ymax
                eosPointer->UpdateState_TPX(props_refine_check[3], ymin_quad + physical_length_quad[1],     forest->m_constZ, xmin_quad + physical_length_quad[0]);        //xmax,ymax
                eosPointer->UpdateState_TPX(props_refine_check[4], ymin_quad + physical_length_quad[1]/2.0, forest->m_constZ, xmin_quad + physical_length_quad[0]/2.0);    //xc,yc

            }
                break;
            case LOOKUPTABLE_FOREST::CONST_NO_VAR_TorHPX: //3D case, the dim must be 3, because it has been checked in the forest construct function
            {
                double dz_physical_quad = physical_length_quad[2] / (num_sample_z - 1.0);
                double z_physical_quad;
                for (int iz = 0; iz < num_sample_z; iz++)
                {
                    z_physical_quad = zmin_quad + dz_physical_quad*iz;
                    for (int iy = 0; iy < num_sample_y; iy++)
                    {
                        y_physical_quad = ymin_quad + dy_physical_quad*iy;
                        for (int ix = 0; ix < num_sample_x; ix++)
                        {
                            x_physical_quad = xmin_quad + dx_physical_quad*ix;
                            regionIndex[iz*num_sample_x*num_sample_y + iy*num_sample_x + ix] = eosPointer->findPhaseRegion_TPX(x_physical_quad, y_physical_quad, z_physical_quad);
                        }
                    }
                }
                data->phaseRegion_cell     = eosPointer->findPhaseRegion_TPX(xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0, zmin_quad + physical_length_quad[2]/2.0);
                // ========================================================
                // calculate properties: four vertices and one midpoint
                eosPointer->UpdateState_TPX(props_refine_check[0], xmin_quad,                               ymin_quad,                                zmin_quad);      //xmin,ymin
                eosPointer->UpdateState_TPX(props_refine_check[1], xmin_quad + physical_length_quad[0],     ymin_quad,                                zmin_quad);      //xmax,ymin
                eosPointer->UpdateState_TPX(props_refine_check[2], xmin_quad,                               ymin_quad + physical_length_quad[1],      zmin_quad);      //xmin,ymax
                eosPointer->UpdateState_TPX(props_refine_check[3], xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      zmin_quad);      //xmax,ymax
                eosPointer->UpdateState_TPX(props_refine_check[4], xmin_quad,                               ymin_quad,                                zmin_quad + physical_length_quad[2]);      //xmin,ymin
                eosPointer->UpdateState_TPX(props_refine_check[5], xmin_quad + physical_length_quad[0],     ymin_quad,                                zmin_quad + physical_length_quad[2]);      //xmax,ymin
                eosPointer->UpdateState_TPX(props_refine_check[6], xmin_quad,                               ymin_quad + physical_length_quad[1],      zmin_quad + physical_length_quad[2]);      //xmin,ymax
                eosPointer->UpdateState_TPX(props_refine_check[7], xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      zmin_quad + physical_length_quad[2]);      //xmax,ymax
                eosPointer->UpdateState_TPX(props_refine_check[8], xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0,  zmin_quad+ physical_length_quad[2]/2.0);      //xc,yc

            }
                break;
            default:
                throw NotImplementedError ("It is impossible! The forest->m_const_which_var is not one of CONST_X_VAR_TorHP, CONST_TorH_VAR_XP, CONST_P_VAR_XTorH and CONST_NO");
                break;
        }

        // ========== 1. refinement check for phase boundary ============
        bool isSame_phaseIndex = true;
        for (int i = 1; i < num_sample_x*num_sample_y*num_sample_z; i++)
        {
            isSame_phaseIndex = (isSame_phaseIndex && (regionIndex[0] == regionIndex[i]));
            if(!isSame_phaseIndex) //if phase indices of all sample points are equal, do not refine; otherwise do refine
            {
                data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_PhaseBoundary;
                data->phaseRegion_cell = MixPhaseRegion;
                break;
            }
        }
        // call MACRO function calculate property refine criterion
        // \todo what is the best criterion for each property ?
        // ========== 2. refinement check for Rho =====================
        if(data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
             CHECK_REFINE_PROP_RMSD(Rho); //!!!! check it later
        }
        if(data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
            CHECK_REFINE_PROP_RMSD(Mu); //!!!! check it later
        }
        // // ========== 2. refinement check for H =====================
        // if(data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
        //     CHECK_REFINE_PROP_RMSD(Rho);
        // }

        delete[] props_refine_check;
        // ============ return refine indicator ===========
        if(quad->level >= forest->m_max_level)return false; //safe check!!! must be at front
        if(data->need_refine)return true;

        return false;
    }

    template <int dim, typename USER_DATA>
    bool RefineFunc_PHX(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, LOOKUPTABLE_FOREST::Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level)
    {
        // cout<<"level: "<<quad->level<<endl;
        if(quad->isHasChildren) return true; //if a quad has children, of course it need refine, but we don't need do anything at here, just return true.

        cxThermal* eosPointer =(cxThermal*)(forest->m_eosPointer); //read only! please DO NOT use this pointer to change any data in the EOS object!!! although it can change the member data.
        // double xv, xl;
        USER_DATA       *data = (USER_DATA *) quad->qData.leaf->user_data;
        data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_NoNeed; //initialize need_refine to be NoNeed
        // const double* physical_length_quad = forest->m_physical_length_quad[quad->level];
        double physical_length_quad[3];
        forest->get_quadrant_physical_length(quad->level, physical_length_quad);
        const int num_sample_x =2; //if want to make more safe check, could increase this number to generate more sample points, for most of cases, the default value 2 is enough.
        const int num_sample_y =2;
        const int num_sample_z = dim==3 ? 2 : 1;
        double dx_physical_quad = physical_length_quad[0] / (num_sample_x - 1.0);
        double dy_physical_quad = physical_length_quad[1] / (num_sample_y - 1.0);
        double x_physical_quad, y_physical_quad;
        ThermodynamicProperties tmp_prop;
        PhaseRegion regionIndex[num_sample_x*num_sample_y*num_sample_z];
        ThermodynamicProperties* props_refine_check = NULL;
        const bool refine_prop = true;
        if(refine_prop)props_refine_check = new ThermodynamicProperties[forest->m_num_children + 1];
        switch (forest->m_const_which_var)
        {
            case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        eosPointer->UpdateState_HPX(tmp_prop, x_physical_quad, y_physical_quad, forest->m_constZ);
                        regionIndex[iy*num_sample_x + ix] = tmp_prop.phase;
                    }
                }
                eosPointer->UpdateState_HPX(tmp_prop, xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0, forest->m_constZ);
                data->phaseRegion_cell     = tmp_prop.phase;
                // ========================================================
                // calculate properties: four vertices and one midpoint
                if(refine_prop)
                {
                    eosPointer->UpdateState_HPX(props_refine_check[0] , xmin_quad,                               ymin_quad,                                forest->m_constZ);      //xmin,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[1] , xmin_quad + physical_length_quad[0],     ymin_quad,                                forest->m_constZ);      //xmax,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[2] , xmin_quad,                               ymin_quad + physical_length_quad[1],      forest->m_constZ);      //xmin,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[3] , xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      forest->m_constZ);      //xmax,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[4] , xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0,  forest->m_constZ);      //xc,yc
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        eosPointer->UpdateState_HPX(tmp_prop, forest->m_constZ, y_physical_quad, x_physical_quad);
                        regionIndex[iy*num_sample_x + ix] = tmp_prop.phase;
                    }
                }
                eosPointer->UpdateState_HPX(tmp_prop, forest->m_constZ, ymin_quad + physical_length_quad[1]/2.0, xmin_quad + physical_length_quad[0]/2.0);
                data->phaseRegion_cell     = tmp_prop.phase;
                // ========================================================
                // calculate properties: four vertices and one midpoint
                if(refine_prop)
                {
                    eosPointer->UpdateState_HPX(props_refine_check[0] , ymin_quad,                                forest->m_constZ,   xmin_quad);                                  //xmin,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[1] , ymin_quad,                                forest->m_constZ,   xmin_quad + physical_length_quad[0]);        //xmax,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[2] , ymin_quad + physical_length_quad[1],      forest->m_constZ,   xmin_quad);                                  //xmin,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[3] , ymin_quad + physical_length_quad[1],      forest->m_constZ,   xmin_quad + physical_length_quad[0]);        //xmax,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[4] , ymin_quad + physical_length_quad[1]/2.0,  forest->m_constZ,   xmin_quad + physical_length_quad[0]/2.0);    //xc,yc
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
            {
                for (int iy = 0; iy < num_sample_y; iy++)
                {
                    y_physical_quad = ymin_quad + dy_physical_quad*iy;
                    for (int ix = 0; ix < num_sample_x; ix++)
                    {
                        x_physical_quad = xmin_quad + dx_physical_quad*ix;
                        eosPointer->UpdateState_HPX(tmp_prop, y_physical_quad, forest->m_constZ, x_physical_quad);
                        regionIndex[iy*num_sample_x + ix] = tmp_prop.phase;
                    }
                }
                eosPointer->UpdateState_HPX(tmp_prop, ymin_quad + physical_length_quad[1]/2.0,  forest->m_constZ, xmin_quad + physical_length_quad[0]/2.0);
                data->phaseRegion_cell     = tmp_prop.phase;
                // ========================================================
                // calculate properties: four vertices and one midpoint
                if(refine_prop)
                {
                    eosPointer->UpdateState_HPX(props_refine_check[0] , ymin_quad,                               forest->m_constZ, xmin_quad);                                  //xmin,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[1] , ymin_quad,                               forest->m_constZ, xmin_quad + physical_length_quad[0]);        //xmax,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[2] , ymin_quad + physical_length_quad[1],     forest->m_constZ, xmin_quad);                                  //xmin,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[3] , ymin_quad + physical_length_quad[1],     forest->m_constZ, xmin_quad + physical_length_quad[0]);        //xmax,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[4] , ymin_quad + physical_length_quad[1]/2.0, forest->m_constZ, xmin_quad + physical_length_quad[0]/2.0);    //xc,yc
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_NO_VAR_TorHPX: //3D case, the dim must be 3, because it has been checked in the forest construct function
            {
                double dz_physical_quad = physical_length_quad[2] / (num_sample_z - 1.0);
                double z_physical_quad;
                for (int iz = 0; iz < num_sample_z; iz++)
                {
                    z_physical_quad = zmin_quad + dz_physical_quad*iz;
                    for (int iy = 0; iy < num_sample_y; iy++)
                    {
                        y_physical_quad = ymin_quad + dy_physical_quad*iy;
                        for (int ix = 0; ix < num_sample_x; ix++)
                        {
                            x_physical_quad = xmin_quad + dx_physical_quad*ix;
                            eosPointer->UpdateState_HPX(tmp_prop, x_physical_quad, y_physical_quad, z_physical_quad);
                            regionIndex[iz*num_sample_x*num_sample_y + iy*num_sample_x + ix] = tmp_prop.phase;
                        }
                    }
                }
                eosPointer->UpdateState_HPX(tmp_prop, xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0, zmin_quad + physical_length_quad[2]/2.0);
                data->phaseRegion_cell     = tmp_prop.phase;
                // ========================================================
                // calculate properties: four vertices and one midpoint
                if(refine_prop)
                {
                    eosPointer->UpdateState_HPX(props_refine_check[0] , xmin_quad,                               ymin_quad,                                zmin_quad);      //xmin,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[1] , xmin_quad + physical_length_quad[0],     ymin_quad,                                zmin_quad);      //xmax,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[2] , xmin_quad,                               ymin_quad + physical_length_quad[1],      zmin_quad);      //xmin,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[3] , xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      zmin_quad);      //xmax,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[4] , xmin_quad,                               ymin_quad,                                zmin_quad + physical_length_quad[2]);      //xmin,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[5] , xmin_quad + physical_length_quad[0],     ymin_quad,                                zmin_quad + physical_length_quad[2]);      //xmax,ymin
                    eosPointer->UpdateState_HPX(props_refine_check[6] , xmin_quad,                               ymin_quad + physical_length_quad[1],      zmin_quad + physical_length_quad[2]);      //xmin,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[7] , xmin_quad + physical_length_quad[0],     ymin_quad + physical_length_quad[1],      zmin_quad + physical_length_quad[2]);      //xmax,ymax
                    eosPointer->UpdateState_HPX(props_refine_check[8] , xmin_quad + physical_length_quad[0]/2.0, ymin_quad + physical_length_quad[1]/2.0,  zmin_quad+ physical_length_quad[2]/2.0);      //xc,yc
                }
            }
                break;
            default:
            throw NotImplementedError("It is impossible! The forest->m_const_which_var is not one of CONST_X_VAR_TorHP, CONST_TorH_VAR_XP, CONST_P_VAR_XTorH and CONST_NO");
                break;
        }

        // ========== 1. refinement check for phase boundary ============
        bool isSame_phaseIndex = true;
        for (int i = 1; i < num_sample_x*num_sample_y*num_sample_z; i++)
        {
            isSame_phaseIndex = (isSame_phaseIndex && (regionIndex[0] == regionIndex[i]));
            if(!isSame_phaseIndex) //if phase indices of all sample points are equal, do not refine; otherwise do refine
            {
                data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_PhaseBoundary;
                data->phaseRegion_cell = MixPhaseRegion;
                break;
            }
        }
        // call MACRO function calculate property refine criterion
        // \todo what is the best criterion for each property ?
        // // ========== 2. refinement check for Rho =====================
        if(refine_prop  && data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
             CHECK_REFINE_PROP_RMSD(Rho);
        }
        if(refine_prop  && data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
            CHECK_REFINE_PROP_RMSD(Mu); //!!!! check it later
        }
        // // ========== 2. refinement check for H =====================
        // if(refine_prop  && data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed){
        //     CHECK_REFINE_PROP_RMSD(H);
        // }


        if(props_refine_check) delete[] props_refine_check;
        // ============ return refine indicator ===========
        if(quad->level >= forest->m_max_level)return false; //safe check!!! must be at front
        if(data->need_refine)return true;

        return false;
    }

    template <int dim, typename USER_DATA>
    void cal_prop_PTX(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, std::map<LOOKUPTABLE_FOREST::Quad_index, LOOKUPTABLE_FOREST::int_pointIndex>& map_ijk2data, double** data)
    {
        cxThermal* eosPointer =(cxThermal*)(forest->m_eosPointer); //read only! please DO NOT use this pointer to change any data in the EOS object!!! although it can change the member data.
        if(eosPointer->get_UpdateWhichProps().size()==0) return; //if no properties are specificed, doesn't do anything and return
        // for parallel computing
        LOOKUPTABLE_FOREST::Quad_index* p_ijk_first = new LOOKUPTABLE_FOREST::Quad_index[map_ijk2data.size()];
        LOOKUPTABLE_FOREST::int_pointIndex* p_ijk_second = new unsigned int[map_ijk2data.size()];
        int num_points=0;
        for(auto &m : map_ijk2data)
        {
            p_ijk_first[num_points] = m.first;
            p_ijk_second[num_points] = m.second;
            num_points++;
        }
        // ----------------------
        double x, y, z;
        ThermodynamicProperties tmp_prop;
        switch (forest->m_const_which_var)
        {
            case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
            {
                // for(auto &ijk2data : map_ijk2data)
                #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    // forest->ijk2xyz(&ijk2data.first, x, y, z);
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_TPX(tmp_prop, x, y, forest->m_constZ);
                    // fill_prop2data(eosPointer, &tmp_prop, eosPointer->m_update_which_props, ijk2data.second);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);

                }
                // std::cout<<"prop calculation done"<<std::endl;
                // int test_id = 0;
                // for(auto &ijk2data : map_ijk2data)
                // {
                //     std::cout<<eosPointer->m_supported_props[eosPointer->m_update_which_props[test_id]]<<ijk2data.second[test_id]<<std::endl;
                // }

            }
                break;
            case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
            {
                //\todo 这里可以考虑并行
                #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_TPX(tmp_prop, forest->m_constZ, y, x);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
            {
            #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_TPX(tmp_prop, y, forest->m_constZ, x);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_NO_VAR_TorHPX: //3D case, the dim must be 3, because it has been checked in the forest construct function
            {
            #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_TPX(tmp_prop, x, y, z);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            default:
                throw NotImplementedError("It is impossible! The forest->m_const_which_var is not one of CONST_X_VAR_TorHP, CONST_TorH_VAR_XP, CONST_P_VAR_XTorH and CONST_NO");
                break;
        }

        // release ijk pointer
        delete[] p_ijk_first;
        delete[] p_ijk_second;
    }

    template <int dim, typename USER_DATA>
    void cal_prop_PHX(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, std::map<LOOKUPTABLE_FOREST::Quad_index, LOOKUPTABLE_FOREST::int_pointIndex>& map_ijk2data, double** data)
    {
        cxThermal* eosPointer =(cxThermal*)(forest->m_eosPointer); //read only! please DO NOT use this pointer to change any data in the EOS object!!! although it can change the member data.
        if(eosPointer->get_UpdateWhichProps().size()==0) return; //if no properties are specificed, doesn't do anything and return
        // for parallel computing
        LOOKUPTABLE_FOREST::Quad_index* p_ijk_first = new LOOKUPTABLE_FOREST::Quad_index[map_ijk2data.size()];
        LOOKUPTABLE_FOREST::int_pointIndex* p_ijk_second = new unsigned int[map_ijk2data.size()];
        int num_points=0;
        for(auto &m : map_ijk2data)
        {
            p_ijk_first[num_points] = m.first;
            p_ijk_second[num_points] = m.second;
            num_points++;
        }
        // ----------------------
        double x, y, z;
        ThermodynamicProperties tmp_prop;
        switch (forest->m_const_which_var)
        {
            case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
            {
                #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_HPX(tmp_prop, x, y, forest->m_constZ);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
                // std::cout<<"prop calculation done"<<std::endl;
                // int test_id = 0;
                // for(auto &ijk2data : map_ijk2data)
                // {
                //     std::cout<<eosPointer->m_supported_props[eosPointer->m_update_which_props[test_id]]<<ijk2data.second[test_id]<<std::endl;
                // }

            }
                break;
            case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
            {
            #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_HPX(tmp_prop, forest->m_constZ, y, x);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
            {
            #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_HPX(tmp_prop, y, forest->m_constZ, x);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_NO_VAR_TorHPX: //3D case, the dim must be 3, because it has been checked in the forest construct function
            {
            #pragma omp parallel for shared(p_ijk_first, p_ijk_second, eosPointer,forest) private(x, y, z, tmp_prop)
                for(int i=0; i<num_points; i++)
                {
                    forest->ijk2xyz(&p_ijk_first[i], x, y, z);
                    eosPointer->UpdateState_HPX(tmp_prop, x, y, z);
                    fill_prop2data(eosPointer, &tmp_prop, eosPointer->get_UpdateWhichProps(), data[p_ijk_second[i]]);
                }
            }
                break;
            default:
                throw NotImplementedError("It is impossible! The forest->m_const_which_var is not one of CONST_X_VAR_TorHP, CONST_TorH_VAR_XP, CONST_P_VAR_XTorH and CONST_NO");
                break;
        }

        // realease p_ijk*
        delete[] p_ijk_first;
        delete[] p_ijk_second;
    }

    /**
     * @brief Uniform refinement. Used to define the minimum level mesh.
     *
     * @tparam dim
     * @tparam USER_DATA
     * @param forest
     * @param quad
     * @param max_level
     * @return true
     * @return false
     */
    template <int dim, typename USER_DATA>
    static bool refine_uniform(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, LOOKUPTABLE_FOREST::Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level)
    {
        if(quad->level < forest->m_min_level)return true;

        return false;
    };

}

#endif