#include "thermo.h"

// ----- LUT related head filess --------
#include "LookUpTableForestI.H"
#include "interpolationI.H"
#include "AMR_LUT_RefineFuncI.H"
//---------------------------------------

namespace xThermal
{
    using namespace std;
    void cxThermal::initialize_data() {
        m_num_threads = 1;
        set_num_threads(m_num_threads);
        m_dim_lut = 0;
        m_pLUT = NULL;
        m_dim_lut_lookup = 0;
        m_pLUT_lookup = NULL;
        init_supported_props();
        m_index_Rho_in_LUT = -1;
        // show progressbar
        m_isShowProgressBar = false;
    }
    cxThermal::cxThermal(/* args */)
    {
        initialize_data();
    }
    
    cxThermal::~cxThermal()
    {
        if(m_pLUT)destroyLUT(m_pLUT, m_dim_lut);
        if(m_pLUT_lookup)destroyLUT(m_pLUT_lookup, m_dim_lut_lookup);
    }

    void cxThermal::set_num_threads(int num_threads)
    {
#if USE_OMP == 1
        omp_set_num_threads(num_threads);
        m_num_threads = num_threads;
#endif
    }
    int cxThermal::get_num_threads()
    {
// #if USE_OMP == 1
//         m_num_threads = omp_get_num_threads();
// #endif
        return  m_num_threads;
    }

//    template<class T> T cxThermal::max_vector(const std::vector<T> &x)
//    {
//        T max = -1E30;
//        std::size_t N = x.size();
//        for (std::size_t i = 0; i < N; ++i)
//        {
//            T axi = x[i];
//            if (axi > max){ max = axi; }
//        }
//        return max;
//    }
//
//    template<class T> T cxThermal::min_vector(const std::vector<T> &x)
//    {
//        T min = 1e40;
//        std::size_t N = x.size();
//        for (std::size_t i = 0; i < N; ++i)
//        {
//            T axi = x[i];
//            if (axi < min){ min = axi; }
//        }
//        return min;
//    }

    /**
     * Triangulate a polygon using the Triangle code, see https://www.cs.cmu.edu/~quake/triangle.html
     * @param x_poly
     * @param y_poly
     * @param trimesh
     */
    void cxThermal::Triangulation(const std::vector<double> &x_poly, const std::vector<double> &y_poly, const double pointInMesh[2], const double dxdy[2], TriMesh &trimesh, bool normalize)
    {
        triangulateio in{}, out{}, out_refine{};

        /* Define input points. */

        in.numberofpoints = (int)x_poly.size();
        in.numberofpointattributes = 1;
        in.numberofsegments = (int)x_poly.size(); //connect
        in.numberofholes = 0;
        in.numberofregions = 1;
        in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
        in.pointattributelist = (REAL *) malloc(in.numberofpoints * in.numberofpointattributes * sizeof(REAL));
        in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int ));
        in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
        // calculate normalize scale
        double scale_x = 1, scale_y = 1;
        if (normalize)
        {
            double x_min = min_vector(x_poly), x_max = max_vector(x_poly);
            double y_min = min_vector(y_poly), y_max = max_vector(y_poly);
            scale_x = 1.0/(x_max - x_min);
            scale_y = 1.0/(y_max - y_min);
        }
        for (int i = 0; i < x_poly.size(); ++i) {
            in.pointlist[0+i*2] = (x_poly[i])*scale_x;
            in.pointlist[1+i*2] = (y_poly[i])*scale_y;
            in.pointattributelist[i] = 0.0;
            in.pointmarkerlist[i] = i;
            in.segmentlist[0+i*2] = i;  //The x_poly,y_poly already have the default connected order!!!
            in.segmentlist[1+i*2] = i+1;
            // std::cout<<"Point: "<<in.pointlist[0+i*2]<<" "<<in.pointlist[1+i*2]<<std::endl;
        }
        in.segmentlist[2*in.numberofsegments-1] = 0; // the last one, connect to the first point
        in.regionlist[0] = (pointInMesh[0])*scale_x;
        in.regionlist[1] = (pointInMesh[1])*scale_y;
        in.regionlist[2] = 1;            /* Regional attribute (for whole mesh). */
        double min_area = dxdy[0]*scale_x * dxdy[1]*scale_y;
        in.regionlist[3] = min_area;          /* Area constraint that will not be used. */

        // printf("Input point set:\n\n");
        // report(&in, 1, 0, 0, 0, 0, 0);
        out.pointlist = (REAL *) nullptr;            /* Not needed if -N switch used. */
        /* Not needed if -N switch used or number of point attributes is zero: */
        out.pointattributelist = (REAL *) nullptr;
        out.pointmarkerlist = (int *) nullptr; /* Not needed if -N or -B switch used. */
        out.trianglelist = (int *) nullptr;          /* Not needed if -E switch used. */
        /* Not needed if -E switch used or number of triangle attributes is zero: */
        out.triangleattributelist = (REAL *) NULL;
        out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
        /* Needed only if segments are output (-p or -c) and -P not used: */
        out.segmentlist = (int *) NULL;
        /* Needed only if segments are output (-p or -c) and -P and -B not used: */
        out.segmentmarkerlist = (int *) NULL;
        out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
        out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
        // triangulate: q20 means the Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.
        triangulate("pQIq20zAena", &in, &out, (struct triangulateio *) NULL);

        //====================refine =================
        /* Needed only if -r and -a switches used: */
        out.trianglearealist = (REAL *) malloc(out.numberoftriangles * sizeof(REAL));
        int ind_vertex = 0;
        for (int i = 0; i < out.numberoftriangles; ++i) {
            for (int j = 0; j < out.numberofcorners; j++) {
                ind_vertex = out.trianglelist[i * out.numberofcorners + j];
                out.trianglearealist[i] = -1;
                // if(out.pointlist[ind_vertex * 2] > -0.1) //10^{xcenter}
                // {
                //     // out.trianglearealist[i] = min_area/10.0;
                //     // continue;
                // }
            }
        }
        out_refine.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
        /* Not needed if -N switch used or number of attributes is zero: */
        out_refine.pointattributelist = (REAL *) NULL;
        out_refine.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
        /* Not needed if -E switch used or number of triangle attributes is zero: */
        out_refine.triangleattributelist = (REAL *) NULL;
        /* Refine the triangulation according to the attached */
        /*   triangle area constraints.                       */
        triangulate("prazBPQ", &out, &out_refine, (struct triangulateio *) NULL);
        //============================================
        //assemble to triMesh
        for (int i = 0; i < out_refine.numberofpoints; i++) {
            trimesh.x.push_back(out_refine.pointlist[i * 2 + 0]/scale_x);
            trimesh.y.push_back(out_refine.pointlist[i * 2 + 1]/scale_y);
            trimesh.z.push_back(0);
        }
        for (int i = 0; i < out_refine.numberoftriangles; i++) {
            std::vector<int> triangle;
            for (int j = 0; j < out_refine.numberofcorners; j++) {
                triangle.push_back(out_refine.trianglelist[i * out_refine.numberofcorners + j]);
            }
            trimesh.connection.push_back(triangle);
        }

        /* Free all allocated arrays, including those allocated by Triangle. */
        if(in.pointlist) free(in.pointlist);
        if(in.pointattributelist) free(in.pointattributelist);
        if(in.pointmarkerlist) free(in.pointmarkerlist);
        if(in.regionlist) free(in.regionlist);
        if(in.segmentlist) free(in.segmentlist);
        if(out.pointlist) free(out.pointlist);
        if(out.pointattributelist) free(out.pointattributelist);
        if(out.pointmarkerlist) free(out.pointmarkerlist);
        if(out.trianglelist) free(out.trianglelist);
        if(out.triangleattributelist) free(out.triangleattributelist);
        // if(out.trianglearealist) free(out.trianglearealist);
        if(out.neighborlist) free(out.neighborlist);
        if(out.segmentlist) free(out.segmentlist);
        if(out.segmentmarkerlist) free(out.segmentmarkerlist);
        if(out.edgelist) free(out.edgelist);
        if(out.edgemarkerlist) free(out.edgemarkerlist);
        if(out_refine.pointlist)free(out_refine.pointlist);
        if(out_refine.pointattributelist)free(out_refine.pointattributelist);
        if(out_refine.trianglelist)free(out_refine.trianglelist);
        if(out_refine.triangleattributelist)free(out_refine.triangleattributelist);
    }

    void cxThermal::writeTriMesh2Txt(const TriMesh &triMesh, std::string path_out) {
        //write to file
        std::ofstream fpout_xy(path_out+"/points.xyz");
        if(!fpout_xy.good()) ERROR("Open file failed: "+path_out+"/points.xyz");
        std::ofstream fpout_connect(path_out+"/connection.txt");
        if(!fpout_connect.good()) ERROR("Open file failed: "+path_out+"/connection.txt");
        for (int i = 0; i < triMesh.x.size(); i++) {
            fpout_xy<<triMesh.x[i]<<" "<<triMesh.y[i]<<" "<<triMesh.z[i]<<std::endl;
        }
        for (int i = 0; i < triMesh.connection.size(); i++) {
            for (int j = 0; j < triMesh.connection[i].size(); j++) {
                fpout_connect<<triMesh.connection[i][j]<<" ";
            }
            fpout_connect<<std::endl;
        }
        fpout_xy.close();
        fpout_connect.close();
        STATUS("The triangle mesh files(points.xyz, connection.txt) have been saved to: "+path_out);
    }
    /**
    * Write a grid mesh to vtu file. The grid is described by three/to 2-D array.
    * @param vtuFile
    * @param XX
    * @param YY
    * @param ZZ
    * @param scale_x
    * @param scale_y
    * @param scale_z
    */
    void cxThermal::writeXXYYZZ2VTU(std::string vtuFile, const std::vector<std::vector<double> > &XX, const std::vector<std::vector<double> > &YY,
                                   const std::vector<std::vector<double> > &ZZ, const double scale_x, const double scale_y,
                                   double scale_z) {
        int nrows = XX.size();
        int ncols = XX[0].size();
        int npoints = nrows * ncols;
        int nCells=(nrows - 1) * (ncols - 1);
        const int VTK_CELLTYPE=9; //quad
        const int np_per_cell=4;
        // write vtu
        std::ofstream fpout(vtuFile);
        if(!fpout.good()) ERROR("Open file failed: "+vtuFile);
        fpout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        fpout<<"  <UnstructuredGrid>\n";
        fpout<<"    <Piece NumberOfPoints=\""<<npoints<<"\" NumberOfCells=\""<<nCells<<"\">\n";
        fpout<<"      <PointData>\n";
        // fpout<<"        <DataArray type=\"Float64\" Name="<<dataname<<" format=\"ascii\">\n";
        // fpout<<"          ";
        // for i in range(0,len(data)):
        //      fpout<<"%f '%(data[i]))
        //  fpout<<"\n        </DataArray>\n";
        fpout<<"      </PointData>\n";
        fpout<<"      <CellData>\n";
        fpout<<"      </CellData>\n";
        fpout<<"      <Points>\n";
        fpout<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j) {
                fpout<<"          "<<XX[i][j]*scale_x<<" "<<YY[i][j]*scale_y<<" "<<ZZ[i][j]*scale_z<<std::endl;
            }
        fpout<<"        </DataArray>\n";
        fpout<<"      </Points>\n";
        fpout<<"      <Cells>\n";
        fpout<<"        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (int nrow = 0; nrow < nrows - 1; ++nrow) {
            for (int ncol = 0; ncol < ncols-1; ++ncol) {
                int LL=ncol + nrow*ncols;
                fpout<<"          "<<LL<<" "<<LL+1<<" "<<LL+1+ncols<<" "<<LL+ncols<<std::endl;
            }
        }
        fpout<<"        </DataArray>\n";
        fpout<<"        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        fpout<<"          ";
        for (int i = 0; i < nCells; ++i) {
            fpout<<(i+1)*np_per_cell<<" ";
        }
        fpout<<"\n        </DataArray>\n";
        fpout<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        fpout<<"          ";
        for (int i = 0; i < nCells; ++i) {
            fpout<<VTK_CELLTYPE<<" ";
        }
        fpout<<"\n        </DataArray>\n";
        fpout<<"      </Cells>\n";
        fpout<<"    </Piece>\n";
        fpout<<"  </UnstructuredGrid>\n";
        fpout<<"</VTKFile>\n";
        fpout.close();
    }

    void cxThermal::writeLine2VTU(std::string vtuFile, const std::vector<double> &X, const std::vector<double> &Y,
                                 const std::vector<double> &Z, const double scale_x, const double scale_y,
                                 double scale_z) {
        int npoints = X.size();
        int nCells = 1; //single line
        const int VTK_CELLTYPE=4; //VTK_POLY_LINE
        // write vtu
        std::ofstream fpout(vtuFile);
        if(!fpout.good()) ERROR("Open file failed: "+vtuFile);
        fpout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        fpout<<"  <UnstructuredGrid>\n";
        fpout<<"    <Piece NumberOfPoints=\""<<npoints<<"\" NumberOfCells=\""<<1<<"\">\n";
        fpout<<"      <PointData>\n";
        // fpout<<"        <DataArray type=\"Float64\" Name="<<dataname<<" format=\"ascii\">\n";
        // fpout<<"          ";
        // for i in range(0,len(data)):
        //      fpout<<"%f '%(data[i]))
        //  fpout<<"\n        </DataArray>\n";
        fpout<<"      </PointData>\n";
        fpout<<"      <CellData>\n";
        fpout<<"      </CellData>\n";
        fpout<<"      <Points>\n";
        fpout<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < npoints; ++i)fpout<<"          "<<X[i]*scale_x<<" "<<Y[i]*scale_y<<" "<<Z[i]*scale_z<<std::endl;
        fpout<<"        </DataArray>\n";
        fpout<<"      </Points>\n";
        fpout<<"      <Cells>\n";
        fpout<<"        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        fpout<<"          ";
        for (int i = 0; i < npoints; ++i)fpout<<i<<" ";
        fpout<<"\n        </DataArray>\n";
        fpout<<"        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        fpout<<"          ";
        fpout<<npoints;
        fpout<<"\n        </DataArray>\n";
        fpout<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        fpout<<"          ";
        for (int i = 0; i < nCells; ++i) {
            fpout<<VTK_CELLTYPE<<" ";
        }
        fpout<<"\n        </DataArray>\n";
        fpout<<"      </Cells>\n";
        fpout<<"    </Piece>\n";
        fpout<<"  </UnstructuredGrid>\n";
        fpout<<"</VTKFile>\n";
        fpout.close();
    }
    void cxThermal::writePhaseBoundaries2VTU(std::string outputPath, const PhaseBoundaries& phaseBoundaries, double scale_T, double scale_p, double scale_X) {
        // save a vtm file
        std::string vtmfile = outputPath+"/phaseBoundary.vtm";
        std::ofstream fpout_vtm(vtmfile);
        if(!fpout_vtm.good()) ERROR("Open file failed: "+vtmfile);
        fpout_vtm<<"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        fpout_vtm<<"  <vtkMultiBlockDataSet>\n";
        // 1. write surface
        int nsurface = phaseBoundaries.surfaces.size();
        for (int i = 0; i < nsurface; ++i) {
            std::string fname = phaseBoundaries.surfaces[i].shortName+".vtu";
            std::string vtufile = outputPath+"/"+fname;
            writeXXYYZZ2VTU(vtufile, phaseBoundaries.surfaces[i].mesh.X, phaseBoundaries.surfaces[i].mesh.T, phaseBoundaries.surfaces[i].mesh.p, scale_X, scale_T, scale_p);
            fpout_vtm<<"    <Block index=\""<<i<<"\" name=\""<<phaseBoundaries.surfaces[i].name<<"\">\n";
            fpout_vtm<<"      <DataSet index=\""<<0<<"\" name=\""<<1<<"\" file=\""<<fname<<"\">\n";
            fpout_vtm<<"      </DataSet>\n";
            fpout_vtm<<"    </Block>\n";
        }
        fpout_vtm<<"  </vtkMultiBlockDataSet>\n";
        fpout_vtm<<"</VTKFile>";
        fpout_vtm.close();
        //2. write lines
        int nlines = phaseBoundaries.lines.size();
        for (int i = 0; i < nlines; ++i) {
            std::string fname = phaseBoundaries.lines[i].shortName+".vtu";
            std::string vtufile = outputPath+"/"+fname;
            writeLine2VTU(vtufile, phaseBoundaries.lines[i].X, phaseBoundaries.lines[i].T, phaseBoundaries.lines[i].p, scale_X, scale_T, scale_p);
        }
        //2. write points to a single vtu file
        int nPoints = phaseBoundaries.points.size();
        std::vector<double> x_points, y_points, z_points;
        for (int i = 0; i < nPoints; ++i) {
            x_points.push_back(phaseBoundaries.points[i].X);
            y_points.push_back(phaseBoundaries.points[i].T);
            z_points.push_back(phaseBoundaries.points[i].p);
        }
        std::string vtufile = outputPath+"/points.vtu";
        writeLine2VTU(vtufile, x_points,y_points,z_points, scale_X, scale_T, scale_p);

    }

    /**
    * Normalize phase boundaries p,T,X to [0,1] range. This function is designed for VTU format output, because the data visualization in ParaView is default supported x,y,z in 1:1:1 scale.
    * @param phaseBoundaries
    */
    void cxThermal::normalizePhaseBoundaries(PhaseBoundaries &phaseBoundaries) {
        double len_p = phaseBoundaries.pmax - phaseBoundaries.pmin;
        double len_T = phaseBoundaries.Tmax - phaseBoundaries.Tmin;
        double len_X = phaseBoundaries.Xmax - phaseBoundaries.Xmin;
        double Xmin = phaseBoundaries.Xmin;
        switch (phaseBoundaries.scale_X) {
            case SCALE_X_log:
            {
                len_X = log10(phaseBoundaries.Xmax) - log10(phaseBoundaries.Xmin);
                Xmin = log10(phaseBoundaries.Xmin);
                // surface
                for (int iSurface = 0; iSurface < phaseBoundaries.surfaces.size(); ++iSurface) {
                    for (int row = 0; row < phaseBoundaries.surfaces[iSurface].mesh.p.size(); ++row) {
                        for (int col = 0; col < phaseBoundaries.surfaces[iSurface].mesh.p[0].size(); ++col) {
                            phaseBoundaries.surfaces[iSurface].mesh.X[row][col] = log10(phaseBoundaries.surfaces[iSurface].mesh.X[row][col]);
                        }
                    }
                }
                // line
                for (int iLine = 0; iLine < phaseBoundaries.lines.size(); ++iLine) {
                    for (int i = 0; i < phaseBoundaries.lines[iLine].p.size(); ++i) {
                        phaseBoundaries.lines[iLine].X[i] = log10(phaseBoundaries.lines[iLine].X[i]);
                    }
                }
                // point
                for (int iPoint = 0; iPoint < phaseBoundaries.points.size(); ++iPoint) {
                    phaseBoundaries.points[iPoint].X = log10(phaseBoundaries.points[iPoint].X);
                }
            }
                break;
            case SCALE_X_loglinear:
            {
                len_X = 1 + phaseBoundaries.ratio_log_to_linear; // the xmax of the map result always = 1, the map range is [-ratio_log_to_linear, 1]
                Xmin = -phaseBoundaries.ratio_log_to_linear;
                // map log and linear part to the same scale
                double Xmin_log = log10(phaseBoundaries.Xmin);
                double Xmax_log = log10(phaseBoundaries.Xcenter);
                double len_log = (Xmax_log - Xmin_log)/phaseBoundaries.ratio_log_to_linear;
                double len_linear = phaseBoundaries.Xmax - phaseBoundaries.Xcenter;
                // surface
                for (int iSurface = 0; iSurface < phaseBoundaries.surfaces.size(); ++iSurface) {
                    for (int row = 0; row < phaseBoundaries.surfaces[iSurface].mesh.p.size(); ++row) {
                        for (int col = 0; col < phaseBoundaries.surfaces[iSurface].mesh.p[0].size(); ++col) {
                            if (phaseBoundaries.surfaces[iSurface].mesh.X[row][col]<=phaseBoundaries.Xcenter)
                            phaseBoundaries.surfaces[iSurface].mesh.X[row][col] = (log10(phaseBoundaries.surfaces[iSurface].mesh.X[row][col]) - Xmax_log)/len_log;
                            else
                                phaseBoundaries.surfaces[iSurface].mesh.X[row][col] = (phaseBoundaries.surfaces[iSurface].mesh.X[row][col] - phaseBoundaries.Xcenter)/len_linear;
                        }
                    }
                }
                // line
                for (int iLine = 0; iLine < phaseBoundaries.lines.size(); ++iLine) {
                    for (int i = 0; i < phaseBoundaries.lines[iLine].p.size(); ++i) {
                        if(phaseBoundaries.lines[iLine].X[i]<=phaseBoundaries.Xcenter)
                        phaseBoundaries.lines[iLine].X[i] = (log10(phaseBoundaries.lines[iLine].X[i]) - Xmax_log)/len_log;
                        else phaseBoundaries.lines[iLine].X[i] = (phaseBoundaries.lines[iLine].X[i] - phaseBoundaries.Xcenter)/len_linear;
                    }
                }
                // point
                for (int iPoint = 0; iPoint < phaseBoundaries.points.size(); ++iPoint) {
                    if(phaseBoundaries.points[iPoint].X<=phaseBoundaries.Xcenter)
                    phaseBoundaries.points[iPoint].X = (log10(phaseBoundaries.points[iPoint].X) - Xmax_log)/len_log;
                    else phaseBoundaries.points[iPoint].X = (phaseBoundaries.points[iPoint].X - phaseBoundaries.Xcenter)/len_linear;
                }
            }
                break;
        }
        // surface
        for (int iSurface = 0; iSurface < phaseBoundaries.surfaces.size(); ++iSurface) {
            for (int row = 0; row < phaseBoundaries.surfaces[iSurface].mesh.p.size(); ++row) {
                for (int col = 0; col < phaseBoundaries.surfaces[iSurface].mesh.p[0].size(); ++col) {
                    phaseBoundaries.surfaces[iSurface].mesh.p[row][col] = (phaseBoundaries.surfaces[iSurface].mesh.p[row][col]-phaseBoundaries.pmin)/len_p;
                    phaseBoundaries.surfaces[iSurface].mesh.T[row][col] = (phaseBoundaries.surfaces[iSurface].mesh.T[row][col]-phaseBoundaries.Tmin)/len_T;
                    phaseBoundaries.surfaces[iSurface].mesh.X[row][col] = (phaseBoundaries.surfaces[iSurface].mesh.X[row][col]-Xmin)/len_X;
                }
            }
        }
        // line
        for (int iLine = 0; iLine < phaseBoundaries.lines.size(); ++iLine) {
            for (int i = 0; i < phaseBoundaries.lines[iLine].p.size(); ++i) {
                phaseBoundaries.lines[iLine].p[i] = (phaseBoundaries.lines[iLine].p[i] - phaseBoundaries.pmin)/len_p;
                phaseBoundaries.lines[iLine].T[i] = (phaseBoundaries.lines[iLine].T[i] - phaseBoundaries.Tmin)/len_T;
                phaseBoundaries.lines[iLine].X[i] = (phaseBoundaries.lines[iLine].X[i] - Xmin)/len_X;
            }
        }
        // point
        for (int iPoint = 0; iPoint < phaseBoundaries.points.size(); ++iPoint) {
            phaseBoundaries.points[iPoint].p = (phaseBoundaries.points[iPoint].p - phaseBoundaries.pmin)/len_p;
            phaseBoundaries.points[iPoint].T = (phaseBoundaries.points[iPoint].T - phaseBoundaries.Tmin)/len_T;
            phaseBoundaries.points[iPoint].X = (phaseBoundaries.points[iPoint].X - Xmin)/len_X;
        }
    }

    //LUT related function
    void cxThermal::init_supported_props() {
        // Temperature
        strcpy(m_supported_props[Update_prop_T].longName          , "Temperature");
        strcpy(m_supported_props[Update_prop_T].shortName         , "T");
        strcpy(m_supported_props[Update_prop_T].unit              , "[K]");
        // Bulk density
        strcpy(m_supported_props[Update_prop_Rho].longName          , "Bulk density");
        strcpy(m_supported_props[Update_prop_Rho].shortName         , "Rho");
        strcpy(m_supported_props[Update_prop_Rho].unit              , "[kg/m3]");
        // Specific enhalpy
        strcpy(m_supported_props[Update_prop_H].longName            , "Bulk specific enthalpy");
        strcpy(m_supported_props[Update_prop_H].shortName           , "H");
        strcpy(m_supported_props[Update_prop_H].unit                , "[J/kg]");
        // Specific heat
        strcpy(m_supported_props[Update_prop_Cp].longName            , "Bulk specific heat capacity");
        strcpy(m_supported_props[Update_prop_Cp].shortName           , "Cp");
        strcpy(m_supported_props[Update_prop_Cp].unit                , "[J/kg/K]");
        // Dynamic viscosity
        strcpy(m_supported_props[Update_prop_Mu].longName            , "Dynamic viscosity");
        strcpy(m_supported_props[Update_prop_Mu].shortName           , "Mu");
        strcpy(m_supported_props[Update_prop_Mu].unit                , "[Pa s]");
        // Derivative: Isothermal compressibility
        strcpy(m_supported_props[Update_prop_IsothermalCompressibility].longName       , "Isothermal compressibility");
        strcpy(m_supported_props[Update_prop_IsothermalCompressibility].shortName      , "kappa");
        strcpy(m_supported_props[Update_prop_IsothermalCompressibility].unit           , "[1/Pa]");
        // Isobaric expansivity
        strcpy(m_supported_props[Update_prop_IsobaricExpansivity].longName       , "Isobaric expansivity");
        strcpy(m_supported_props[Update_prop_IsobaricExpansivity].shortName      , "beta");
        strcpy(m_supported_props[Update_prop_IsobaricExpansivity].unit           , "[1/K]");
        // // Derivative: drhodh
        // strcpy(m_supported_props[Update_prop_drhodh].longName       , "dRho/dH");
        // strcpy(m_supported_props[Update_prop_drhodh].shortName      , "dRhodH");
        // strcpy(m_supported_props[Update_prop_drhodh].unit           , "[kg2/(m3 J)]");
    }

    void cxThermal::parse_update_which_props(int update_which_props) {
        if(!m_update_which_props.empty())m_update_which_props.clear(); //safety check

        // bitmask
        for(auto &ind2name_prop : m_supported_props)
        {
            if( (update_which_props & ind2name_prop.first) == ind2name_prop.first)
            {
                m_update_which_props[ind2name_prop.first] = ind2name_prop.second;
            }
        }
        // print info
        STATUS("Update properties: "+to_string(m_update_which_props.size()));
        int ind = 0;
        for (auto &m : m_update_which_props)
        {
            ind++;
            STATUS(to_string(ind) + " : " + m.second.longName);
        }
    }

    void cxThermal::destroyLUT(void* pLUT, int& dim_lut) {
        if(pLUT)
        {
            if(dim_lut==2)
            {
                // static_cast<LookUpTableForest_2D*>(m_lut)->destory();
                // ((LookUpTableForest_2D *)m_pLUT)->destory();
                delete ((LOOKUPTABLE_FOREST::LookUpTableForest_2D *)pLUT);
            }else
            {
                // ((LookUpTableForest_3D *)m_pLUT)->destory();
                delete ((LOOKUPTABLE_FOREST::LookUpTableForest_3D *)pLUT);
            }
            pLUT = NULL;
            dim_lut = 0;
        }
    }

    void cxThermal::save_lut_to_vtk(string filename, bool isNormalizeXYZ) {
        if(m_pLUT)
        {
            if(m_dim_lut==2)((LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT)->write_to_vtk(filename, isNormalizeXYZ);
            else ((LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT)->write_to_vtk(filename, isNormalizeXYZ);
        }else if(m_pLUT_lookup) //load_LUT from file, the pointer is m_pLUT_lookup
        {
            if(m_dim_lut_lookup==2)((LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup)->write_to_vtk(filename, isNormalizeXYZ);
            else ((LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT_lookup)->write_to_vtk(filename, isNormalizeXYZ);
        }
    }

    void cxThermal::save_lut_to_binary(string filename) {
        if(m_pLUT)
        {
            if(m_dim_lut==2)((LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT)->write_to_binary(filename);
            else ((LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT)->write_to_binary(filename);
        }
    }

    void cxThermal::createLUT_2D(double *xy_min, double *xy_max, double constZ,
                               LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var,
                               LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level, int max_level,
                               int update_which_props) {
        // parsing which properties need to be calculated
        parse_update_which_props(update_which_props);

        destroyLUT(m_pLUT, m_dim_lut); //destroy lut pointer and release all data before create a new one.
        // WAIT("destroyLUT");
        clock_t start = clock();
        STATUS("Creating 2D lookup table ...");
        // const int dim =2;
        m_dim_lut = 2;
        LOOKUPTABLE_FOREST::LookUpTableForest_2D* tmp_lut_2D = new LOOKUPTABLE_FOREST::LookUpTableForest_2D (xy_min, xy_max, constZ, const_which_var, TorH, max_level, m_update_which_props, this);
        m_pLUT = tmp_lut_2D;
        // WAIT("new LookUpTableForest_2D");
        // refine
        tmp_lut_2D->set_min_level(min_level);
        tmp_lut_2D->refine(refine_uniform);
        // WAIT("refine_uniform");
        // parallel refine
        if (tmp_lut_2D->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
        {
#ifdef USE_OMP
#pragma omp parallel //shared(n)
#endif
            {
#ifdef USE_OMP
#pragma omp single
#endif
                {
                    printf("Do refinement using %d threads.\n", m_num_threads);
                    tmp_lut_2D->refine(RefineFunc_PTX);
                }
            }
            // update properties data on leaves
            STATUS_time("Lookup table refinement done", (clock() - start)/m_num_threads);
            tmp_lut_2D->construct_props_leaves(cal_prop_PTX);
        }else if(tmp_lut_2D->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
        {
#ifdef USE_OMP
#pragma omp parallel //shared(n)
#endif
            {
#ifdef USE_OMP
#pragma omp single
#endif
                {
                    printf("Do refinement using %d threads.\n", m_num_threads);
                    tmp_lut_2D->refine(RefineFunc_PHX);
                }
            }

            //
            STATUS_time("Lookup table refinement done", (clock() - start)/m_num_threads);
            tmp_lut_2D->construct_props_leaves(cal_prop_PHX);
        }else
        {
            ERROR("The EOS space only support TPX and HPX!");
        }
        tmp_lut_2D->print_summary();
        // WAIT("createLUT_2D");
        // tmp_lut_2D->m_map_ijk2data.clear();
        // WAIT("delete map_ijk2data");
    }

    void cxThermal::createLUT_2D(double xmin, double xmax, double ymin, double ymax, double constZ,
                               LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var,
                               LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level, int max_level,
                               int update_which_props) {
        double xy_min[2] = {xmin, ymin};
        double xy_max[2] = {xmax, ymax};
        createLUT_2D(xy_min, xy_max, constZ, const_which_var, TorH,  min_level, max_level, update_which_props);
    }

    const std::map<int, propInfo> &cxThermal::get_UpdateWhichProps() {
        return m_update_which_props;;
    }
    std::string cxThermal::phase_name(PhaseRegion phase_index)
    {
        return map_phase2name[phase_index];
    }

    ThermodynamicProperties cxThermal::Boiling_p_props(const double &T) {
        ThermodynamicProperties props;
        Boiling_p(T, props);
        return props;
    }

    ThermodynamicProperties cxThermal::Boiling_T_props(const double& p)
    {
        ThermodynamicProperties props;
        Boiling_T(p, props);
        return props;
    }

    LOOKUPTABLE_FOREST::Quadrant<2, LOOKUPTABLE_FOREST::FIELD_DATA<2> > *
    cxThermal::lookup(ThermodynamicProperties &prop, double x, double y) {
        LOOKUPTABLE_FOREST::LookUpTableForest_2D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup; // make temporary copy of the pointer
        // cout<<"constZ: "<<tmp_lut->m_constZ<<", y: "<<y<<", x: "<<x<<endl;
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *targetLeaf = NULL;
        double xyz_min_target[2];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, tmp_lut->m_constZ);
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
            {
                switch (tmp_lut->m_const_which_var)
                {
                    case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                        UpdateState_TPX(prop,x, y, tmp_lut->m_constZ);
                        break;
                    case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                        UpdateState_TPX(prop,y, tmp_lut->m_constZ, x);
                        break;
                    case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                        UpdateState_TPX(prop,tmp_lut->m_constZ, y, x);
                        break;
                    default:
                    ERROR("Impossible case occurs in LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > * cH2ONaCl::lookup(ThermodynamicProperties& prop, double x, double y)");
                        break;
                }
            }else if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
            {
                switch (tmp_lut->m_const_which_var)
                {
                    case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                         UpdateState_HPX(prop,x, y, tmp_lut->m_constZ);
                        break;
                    case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                        UpdateState_HPX(prop,y, tmp_lut->m_constZ, x);
                        break;
                    case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                        UpdateState_HPX(prop,tmp_lut->m_constZ, y, x);
                        break;
                    default:
                    ERROR("Impossible case occurs in LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > * cH2ONaCl::lookup(ThermodynamicProperties& prop, double x, double y)");
                        break;
                }
            }else
            {
                ERROR("The EOS space only support TPX and HPX!");
            }
        }
        else
        {
            double xy[2] = {x, y};
            interp_quad_prop<2>(targetLeaf,xyz_min_target, prop, xy);
        }
        return targetLeaf;
    }

    template<int dim>
    void cxThermal::interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim>> *targetLeaf,
                                    double *xyz_min_target, ThermodynamicProperties &prop, const double xyz[dim]) {
        LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> >* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> >*)m_pLUT_lookup;
        double physical_length[dim]; //physical length of the quad
        double coeff[dim][2];
        const int num_children = tmp_lut->m_num_children;
        // double* values_at_vertices = new double[num_children];
        double values_at_vertices[tmp_lut->m_num_children];
        tmp_lut->get_quadrant_physical_length(targetLeaf->level, physical_length);
        get_coeff_bilinear<dim> (xyz_min_target, physical_length, xyz, coeff);
        // Rho
        for (int i = 0; i < tmp_lut->m_num_children; i++){
            values_at_vertices[i] = tmp_lut->m_props_unique_points_leaves.data[targetLeaf->qData.leaf->index_props[i]][tmp_lut->m_map_prop2index[Update_prop_Rho]];
        }
        bilinear_cal<dim>(coeff, values_at_vertices, prop.Rho);
        // H
        for (int i = 0; i < tmp_lut->m_num_children; i++){
            values_at_vertices[i] = tmp_lut->m_props_unique_points_leaves.data[targetLeaf->qData.leaf->index_props[i]][tmp_lut->m_map_prop2index[Update_prop_H]];
        }
        bilinear_cal<dim>(coeff, values_at_vertices, prop.H);

        // phase region
        prop.phase = targetLeaf->qData.leaf->user_data->phaseRegion_cell;

        // delete[] values_at_vertices;
    }

    template<int dim>
    void cxThermal::interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim>> *targetLeaf,
                                    double *xyz_min_target, double *props, const double xyz[dim]) {
        LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> >* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> >*)m_pLUT_lookup;
        double physical_length[dim]; //physical length of the quad
        double coeff[dim][2];
        const int num_children = tmp_lut->m_num_children;
        double* values_at_vertices = new double[tmp_lut->m_num_node_per_quad];
        double** pNodeData = new double*[tmp_lut->m_num_node_per_quad];
        // LOOKUPTABLE_FOREST::Quad_index *ijk_nodes_quad = new LOOKUPTABLE_FOREST::Quad_index[tmp_lut->m_num_node_per_quad]; // \todo 如果使用二阶插值，则需要更多节点，需要通过cellType进行判断：比如二维情况九点quad，那么需要限制max_level必须小于MAX_FOREST_LEVEL-2，不过这个好办，在构造函数里面判断一下进行安全检查就行
        // tmp_lut->get_ijk_nodes_quadrant(targetLeaf, &targetLeaf->qData.leaf->coord.ijk, tmp_lut->m_num_node_per_quad, ijk_nodes_quad);

        tmp_lut->get_quadrant_physical_length((int)targetLeaf->level, physical_length);
        get_coeff_bilinear<dim> (xyz_min_target, physical_length, xyz, coeff);
        // for (int i = 0; i < dim; i++)
        // {
        //     cout<<"    "<<coeff[i][0]<<" "<<coeff[i][1]<<endl;
        // }

        // interpolate props
        int ind_prop =0;
        for (int i = 0; i < tmp_lut->m_num_node_per_quad; i++){
            pNodeData[i] = tmp_lut->m_props_unique_points_leaves.data[targetLeaf->qData.leaf->index_props[i]];
        }
        for(auto &map_props : tmp_lut->m_map_props)
        {
            // if(ind_prop==0)cout<<" prop: "<<map_props.second.longName<<endl;
            for (int i = 0; i < tmp_lut->m_num_node_per_quad; i++){
                values_at_vertices[i] = pNodeData[i][ind_prop]; //\todo need to be optimized, this map index takes more time
                // if(ind_prop==0)cout<<"    node "<<i<<" : "<<values_at_vertices[i]<<endl;
            }
            bilinear_cal<dim>(coeff, values_at_vertices, props[ind_prop]);
            ind_prop++;
        }
        delete[] values_at_vertices;
        delete[] pNodeData;
    }

    LOOKUPTABLE_FOREST::Quadrant<2, LOOKUPTABLE_FOREST::FIELD_DATA<2> > *
    cxThermal::lookup(double *props, double *xyz_min_target, double x, double y, bool is_cal) {
        LOOKUPTABLE_FOREST::LookUpTableForest_2D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup; // make temporary copy of the pointer
        // safety check: bound check
        if(x<tmp_lut->m_xyz_min[0] || x>tmp_lut->m_xyz_max[0] || y<tmp_lut->m_xyz_min[1] || y>tmp_lut->m_xyz_max[1])
        {
            throw OutOfRangeError ("The lookup point: ("+to_string(x)+", "+to_string(y)+") out of lookup table xy range. T ["
            +std::to_string(tmp_lut->m_xyz_min[0]) + ", "+std::to_string(tmp_lut->m_xyz_max[0])+"], p ["
            +std::to_string(tmp_lut->m_xyz_min[1])+", "+std::to_string(tmp_lut->m_xyz_max[1])+"]");
        }
        // -------------------------------
        // cout<<"constZ: "<<tmp_lut->m_constZ<<", y: "<<y<<", x: "<<x<<endl;
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *targetLeaf = nullptr;
        // double xyz_min_target[2];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, tmp_lut->m_constZ);
        // cout<<xyz_min_target[0]-targetLeaf->qData.leaf->coord.xyz[0]<<"   "<<xyz_min_target[1] - targetLeaf->qData.leaf->coord.xyz[1]<<endl;
        // cout<<"  level: "<<targetLeaf->level<<", x: "<<x<<", y: "<<y<<", quad x: "<<targetLeaf->qData.leaf->coord.xyz[0]<<", quad y: "<<targetLeaf->qData.leaf->coord.xyz[1]<<endl;
        ThermodynamicProperties tmp_prop;
        cout<<"targetLeaf->qData.leaf->user_data->need_refine: "<<targetLeaf->qData.leaf->user_data->need_refine<<endl;
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            if(is_cal) // if set is_cal true means cal properties using acurate equation, otherwise interpolate anyway
            {
                if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
                {
                    switch (tmp_lut->m_const_which_var)
                    {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            UpdateState_TPX(tmp_prop, x, y, tmp_lut->m_constZ);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            // cout<<"cal: "<<x<<", "<<y<<", "<<tmp_lut->m_constZ<<": "<<tmp_prop.Rho<<", "<<props[0]<<endl;
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            UpdateState_TPX(tmp_prop, y, tmp_lut->m_constZ, x);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            UpdateState_TPX(tmp_prop, tmp_lut->m_constZ, y, x);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            break;
                        default:
                        ERROR("Impossible case occurs in LOOKUPTABLE_FOREST::Quadrant<2,H2ONaCl::FIELD_DATA<2> > * cH2ONaCl::lookup(double* props, double* xyz_min_target,  double x, double y)");
                            break;
                    }
                }else if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
                {
                    switch (tmp_lut->m_const_which_var)
                    {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            UpdateState_HPX(tmp_prop, x, y, tmp_lut->m_constZ);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            UpdateState_HPX(tmp_prop, y, tmp_lut->m_constZ, x);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            UpdateState_HPX(tmp_prop, tmp_lut->m_constZ, y, x);
                            fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                            break;
                        default:
                        ERROR("Impossible case occurs in LOOKUPTABLE_FOREST::Quadrant<2,H2ONaCl::FIELD_DATA<2> > * cH2ONaCl::lookup(H2ONaCl::PROP_H2ONaCl& prop, double x, double y)");
                            break;
                    }
                }else
                {
                    ERROR("The EOS space only support TPX and HPX!");
                }
            }else
            {
                double xy[2] = {x, y};
                interp_quad_prop<2>(targetLeaf,xyz_min_target, props, xy);
            }
        }
        else
        {
            double xy[2] = {x, y};
            interp_quad_prop<2>(targetLeaf,xyz_min_target, props, xy);
            // cout<<"lookup: "<<x<<", "<<y<<", rho: "<<props[0]<<endl;
        }
        return targetLeaf;
    }

    LOOKUPTABLE_FOREST::Quadrant<3, LOOKUPTABLE_FOREST::FIELD_DATA<3> > *
    cxThermal::lookup(ThermodynamicProperties &prop, double x, double y, double z) {
        if(m_dim_lut_lookup!=3)ERROR("The dim of the LUT is not 3, but you call the 3D lookup function");

        LOOKUPTABLE_FOREST::LookUpTableForest_3D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT_lookup; // make temporary copy of the pointer
        // safety check: bound check
        if(x<tmp_lut->m_xyz_min[0] || x>tmp_lut->m_xyz_max[0] || y<tmp_lut->m_xyz_min[1] || y>tmp_lut->m_xyz_max[1] || z<tmp_lut->m_xyz_min[2] || z>tmp_lut->m_xyz_max[2])
        {
            ERROR("The lookup point: ("+to_string(x)+", "+to_string(y)+", "+to_string(z)+") out of lookup table xyz range.");
        }
        // -------------------------------
        // cout<<"T: "<<x<<", P: "<<y<<", X: "<<z<<endl;
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
        double xyz_min_target[3];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, z);
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
            {
                UpdateState_TPX(prop,x, y, z); //For 3D case, the order of x,y,z MUST BE TorH, p, X.
            }else if (tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
            {
                UpdateState_HPX(prop, x, y, z); //For 3D case, the order of x,y,z MUST BE TorH, p, X.
            }else
            {
                ERROR("The EOS space only support TPX and HPX!");
            }
        }
        else
        {
            double xyz[3] = {x, y, z};
            interp_quad_prop<3>(targetLeaf,xyz_min_target, prop, xyz);
        }
        return targetLeaf;
    }

    LOOKUPTABLE_FOREST::Quadrant<3, LOOKUPTABLE_FOREST::FIELD_DATA<3> > *
    cxThermal::lookup(double *props, double *xyz_min_target, double x, double y, double z, bool is_cal) {
        if(m_dim_lut_lookup!=3)ERROR("The dim of the LUT is not 3, but you call the 3D lookup function");

        LOOKUPTABLE_FOREST::LookUpTableForest_3D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT_lookup; // make temporary copy of the pointer
        // safety check: bound check
        if(x<tmp_lut->m_xyz_min[0] || x>tmp_lut->m_xyz_max[0] || y<tmp_lut->m_xyz_min[1] || y>tmp_lut->m_xyz_max[1] || z<tmp_lut->m_xyz_min[2] || z>tmp_lut->m_xyz_max[2])
        {
            throw OutOfRangeError("The lookup point: ("+to_string(x)+", "+to_string(y)+", "+to_string(z)+") out of lookup table xyz range.");
        }
        // -------------------------------
        // cout<<"T: "<<x<<", P: "<<y<<", X: "<<z<<endl;
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
        // double xyz_min_target[3];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, z);
        ThermodynamicProperties tmp_prop;
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            if(is_cal)
            {
                if(tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
                {
                     UpdateState_TPX(tmp_prop, x, y, z); //For 3D case, the order of x,y,z MUST BE TorH, p, X.
                    fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                }else if (tmp_lut->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
                {
                    UpdateState_HPX(tmp_prop, x, y, z); //For 3D case, the order of x,y,z MUST BE TorH, p, X.
                    fill_prop2data(this, &tmp_prop, tmp_lut->m_map_props, props);
                }else
                {
                    throw ValueError("The EOS space only support TPX and HPX! tmp_lut->m_TorH: "+std::to_string(tmp_lut->m_TorH));
                }
            }else
            {
                double xyz[3] = {x, y, z};
                interp_quad_prop<3>(targetLeaf,xyz_min_target, props, xyz);
            }
        }
        else
        {
            double xyz[3] = {x, y, z};
            interp_quad_prop<3>(targetLeaf,xyz_min_target, props, xyz);
        }
        return targetLeaf;
    }

    LOOKUPTABLE_FOREST::Quadrant<2, LOOKUPTABLE_FOREST::FIELD_DATA<2> > *
    cxThermal::lookup_only(ThermodynamicProperties &prop, double x, double y) {
        LOOKUPTABLE_FOREST::LookUpTableForest_2D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup; // make temporary copy of the pointer
        // cout<<"constZ: "<<tmp_lut->m_constZ<<", y: "<<y<<", x: "<<x<<endl;
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *targetLeaf = NULL;
        double xyz_min_target[2];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, tmp_lut->m_constZ);
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            //only lookup, so don't do anything if the lookup point in the needRefine-cell
            prop.phase = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
        }
        else
        {
            double xy[2] = {x, y};
            interp_quad_prop<2>(targetLeaf,xyz_min_target, prop, xy);
            prop.phase = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
        }
        return targetLeaf;
    }

    LOOKUPTABLE_FOREST::Quadrant<3, LOOKUPTABLE_FOREST::FIELD_DATA<3> > *
    cxThermal::lookup_only(ThermodynamicProperties &prop, double x, double y, double z) {
        if(m_dim_lut_lookup!=3)ERROR("The dim of the LUT is not 3, but you call the 3D lookup function");

        LOOKUPTABLE_FOREST::LookUpTableForest_3D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT_lookup; // make temporary copy of the pointer
        // cout<<"T: "<<x<<", P: "<<y<<", X: "<<z<<endl;
        LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
        double xyz_min_target[3];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, z);
        if(targetLeaf->qData.leaf->user_data->need_refine)
        {
            //this is a looup-only version, don't do anything if the lookup point in the needRefine-cell
            prop.phase = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
        }
        else
        {
            double xyz[3] = {x, y, z};
            interp_quad_prop<3>(targetLeaf,xyz_min_target, prop, xyz);
            prop.phase = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
        }
        return targetLeaf;
    }

    void cxThermal::loadLUT(const string& filename, bool printStatus) {
        destroyLUT(m_pLUT_lookup, m_dim_lut_lookup); //destroy LUT if it already exists.
        m_dim_lut_lookup = LOOKUPTABLE_FOREST::get_dim_from_binary(filename);

        switch (m_dim_lut_lookup)
        {
            case 2:
            {
                m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_2D(filename, this, printStatus));
                auto* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup;
                int ind_rho = 0;
                for(auto & m : tmp_lut->m_map_props)
                {
                    if( (Update_prop_Rho & m.first) == Update_prop_Rho)
                    {
                        m_index_Rho_in_LUT = ind_rho; //it will be used in Rho_lookup function
                    }
                    ind_rho++;
                }
            }
                break;
            case 3:
            {
                m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_3D(filename, this, printStatus));
                int ind_rho = 0;
                auto* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pLUT_lookup;
                for(auto & m : tmp_lut->m_map_props)
                {
                    if( (Update_prop_Rho & m.first) == Update_prop_Rho)
                    {
                        m_index_Rho_in_LUT = ind_rho; //it will be used in Rho_lookup function
                    }
                    ind_rho++;
                }
            }
                break;
            default:
            throw UnableToLoadError("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+filename);
                break;
        }
        if (m_index_Rho_in_LUT<0)
        {
            throw NotImplementedError("The loaded LUT doesn't contains Density, this is a required property.");
        } else
        {
            if(printStatus)STATUS("Find the valid index of density is "+std::to_string(m_index_Rho_in_LUT)+" in the loaded LUT");
        }
    }

    /**
    * @brief Only lookup density for given T,P,X from loaded LUT. This information can be a good guess of Rho for IAPWS95 solution of Rho.
    * @param Rho_estimate Interpolated density in the quad.
    * @param Rho_min Minimum rho in the target (searched) quad
    * @param Rho_max Maximum rho in the target (searched) quad
    * @param T [K]
    * @param P [Pa]
    */
    PhaseRegion cxThermal::Rho_lookup(double &Rho_estimate, double &Rho_min, double &Rho_max, const double &x, const double &y) {
        static bool LUT_loaded=true;
        if (!m_pLUT_lookup)
        {
            //throw ValueError("Look up table is not loaded, please use loadLUT function load a valid LUT first.");
            if(LUT_loaded) WARNING("Look up table is not loaded, can not get estimated Rho from LUT, please use loadLUT function load a valid LUT first to use this speed up feature. This warning will only display once. ")
            LUT_loaded = false;
            return Unknown;
        }
        LOOKUPTABLE_FOREST::LookUpTableForest_2D* tmp_lut = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pLUT_lookup; // make temporary copy of the pointer
        // cout<<"constZ: "<<tmp_lut->m_constZ<<", y: "<<y<<", x: "<<x<<endl;
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *targetLeaf = NULL;
        double xyz_min_target[2];
        tmp_lut->searchQuadrant(targetLeaf, xyz_min_target, x, y, tmp_lut->m_constZ);
        double xy[2] = {x, y};
        // calculate mean value of Rho
        const int dim = 2;
        double physical_length[dim]; //physical length of the quad
        double** pNodeData = new double*[tmp_lut->m_num_node_per_quad];
        tmp_lut->get_quadrant_physical_length((int)targetLeaf->level, physical_length);
        // interpolate props
        int ind_prop =0;
        for (int i = 0; i < tmp_lut->m_num_node_per_quad; i++){
            pNodeData[i] = tmp_lut->m_props_unique_points_leaves.data[targetLeaf->qData.leaf->index_props[i]];
        }
        Rho_estimate = 0, Rho_min = 1E20, Rho_max = -1E20;
        for (int i = 0; i < tmp_lut->m_num_node_per_quad; i++)
        {
            Rho_estimate += pNodeData[i][m_index_Rho_in_LUT];
            Rho_min = std::min(Rho_min, pNodeData[i][m_index_Rho_in_LUT]);
            Rho_max = std::max(Rho_max, pNodeData[i][m_index_Rho_in_LUT]);
        }
        Rho_estimate = Rho_estimate/tmp_lut->m_num_node_per_quad;

        delete[] pNodeData;
        return targetLeaf->qData.leaf->user_data->phaseRegion_cell;
//        return targetLeaf;
    }

    ThermodynamicProperties cxThermal::UpdateState_TPX(const double &T, const double &p, const double &X) {
        ThermodynamicProperties props;
        UpdateState_TPX(props, T, p, X);
        return props;
    }
    ThermodynamicProperties cxThermal::UpdateState_HPX(const double &H, const double &p, const double &X) {
        ThermodynamicProperties props;
        UpdateState_HPX(props, H, p, X);
        return props;
    }

    // not a member function
    void fill_prop2data(cxThermal* pEOS, const ThermodynamicProperties* prop, const std::map<int, propInfo>& update_which_props, double* data)
    {
        // for (size_t i = 0; i < update_which_props.size(); i++)
        int i = 0;
        for(auto &m : update_which_props)
        {
            switch (m.first)
            {
                case Update_prop_Rho:
                    data[i] = prop->Rho;
                    break;
                case Update_prop_H:
                    data[i] = prop->H;
                    break;
                case Update_prop_T:
                    data[i] = prop->T;
                    break;
                case Update_prop_Cp:
                    data[i] = prop->Cp;
                    break;
                case Update_prop_Mu:
                    data[i] = prop->Mu;
                    break;
                case Update_prop_IsobaricExpansivity:
                    data[i] = prop->IsobaricExpansivity;
                    break;
                case Update_prop_IsothermalCompressibility:
                    data[i] = prop->IsothermalCompressibility;
                    break;
                default:
                WARNING("Unsupported property update: " + string(m.second.longName));
                    data[i] = 0;
                    break;
            }
            i++;
        }
    }

    /**
     *
     * @param file_lut The .bin file.
     */
    Head_AMR_LUT cxThermal::getLutInfo(string file_lut, bool printSummary)
    {
        Head_AMR_LUT headinfo;
        int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(file_lut);
        headinfo.dim = m_dim_lut;
        void* m_pLUT_lookup = NULL;
        switch (m_dim_lut)
        {
            case 2:
            {
                LOOKUPTABLE_FOREST::LookUpTableForest_2D* m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_2D(file_lut, NULL,false));
                if(printSummary)m_pLUT_lookup->print_summary();

                int factor = 1<<m_pLUT_lookup->m_max_level;
                if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
                {
                    headinfo.space = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
                    headinfo.spaceName = "TPX";
                }
                else if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
                {
                    headinfo.space = LOOKUPTABLE_FOREST::EOS_ENERGY_H;
                    headinfo.spaceName = "HPX";
                }
                switch (m_pLUT_lookup->m_const_which_var)
                {
                    headinfo.constWhich = m_pLUT_lookup->m_const_which_var;
                    case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                        {
                            headinfo.constWhich_Name = "P";
                            headinfo.constValue = m_pLUT_lookup->m_constZ; //Pa
                            headinfo.var_names.push_back("X");
                            for (int ii = 0; ii < 2; ++ii) {
                                headinfo.var_ranges.push_back({m_pLUT_lookup->m_xyz_min[ii],m_pLUT_lookup->m_xyz_max[ii]});
                                headinfo.var_maxResolution.push_back((m_pLUT_lookup->m_xyz_max[ii] - m_pLUT_lookup->m_xyz_min[ii])/factor);
                            }
                            if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)headinfo.var_names.push_back("T");
                            else if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)headinfo.var_names.push_back("H");
                        }
                        break;
                    case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                    {
                        if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)headinfo.constWhich_Name = "T";
                        else if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)headinfo.constWhich_Name = "H";
                        headinfo.constValue = m_pLUT_lookup->m_constZ; //K
                        headinfo.var_names.push_back("X");
                        headinfo.var_names.push_back("P");
                        for (int ii = 0; ii < 2; ++ii) {
                            headinfo.var_ranges.push_back({m_pLUT_lookup->m_xyz_min[ii],m_pLUT_lookup->m_xyz_max[ii]});
                            headinfo.var_maxResolution.push_back((m_pLUT_lookup->m_xyz_max[ii] - m_pLUT_lookup->m_xyz_min[ii])/factor);
                        }
                    }
                        break;
                    case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                    {
                        headinfo.constWhich_Name = "X";
                        headinfo.constValue = m_pLUT_lookup->m_constZ; //Pa
                        if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)headinfo.var_names.push_back("T");
                        else if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)headinfo.var_names.push_back("H");
                        for (int ii = 0; ii < 2; ++ii) {
                            headinfo.var_ranges.push_back({m_pLUT_lookup->m_xyz_min[ii],m_pLUT_lookup->m_xyz_max[ii]});
                            headinfo.var_maxResolution.push_back((m_pLUT_lookup->m_xyz_max[ii] - m_pLUT_lookup->m_xyz_min[ii])/factor);
                        }
                        headinfo.var_names.push_back("X");
                    }
                        break;
                    case LOOKUPTABLE_FOREST::CONST_NO_VAR_TorHPX:
                    {
                        headinfo.constWhich_Name = "None";
                        if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)headinfo.var_names.push_back("T");
                        else if(m_pLUT_lookup->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)headinfo.var_names.push_back("H");
                        headinfo.var_names.push_back("P");
                        headinfo.var_names.push_back("X");
                        for (int ii = 0; ii < 3; ++ii) {
                            headinfo.var_ranges.push_back({m_pLUT_lookup->m_xyz_min[ii],m_pLUT_lookup->m_xyz_max[ii]});
                            headinfo.var_maxResolution.push_back((m_pLUT_lookup->m_xyz_max[ii] - m_pLUT_lookup->m_xyz_min[ii])/factor);
                        }
                    }
                        break;
                    default:
                        break;
                }
                headinfo.min_level = m_pLUT_lookup->m_min_level;
                headinfo.max_level = m_pLUT_lookup->m_max_level;
                headinfo.num_leaves = m_pLUT_lookup->get_num_leaves();
                headinfo.num_nodes = m_pLUT_lookup->m_props_unique_points_leaves.num_points;
                headinfo.num_need_refine = m_pLUT_lookup->get_num_need_refine();
                headinfo.num_leaves_nextRefine = headinfo.num_need_refine*(m_pLUT_lookup->m_num_children - 1);
                headinfo.num_props = m_pLUT_lookup->m_props_unique_points_leaves.num_props;
                // props info
                int ind =0;
                for(auto & m : m_pLUT_lookup->m_map_props)
                {
                    headinfo.prop_names.push_back(std::string(m.second.longName) + ":" + std::string(m.second.shortName) + ":" + std::string(m.second.unit));
                    ind++;
                }
                const int dim = 2;
                double byte_forest_leaves = sizeof(LOOKUPTABLE_FOREST::LeafQuad<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim>>) * m_pLUT_lookup->get_num_leaves();
                double byte_forest_nonleaves = sizeof(LOOKUPTABLE_FOREST::NonLeafQuad<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim>>) * (m_pLUT_lookup->get_num_quads() - m_pLUT_lookup->get_num_leaves());
                double byte_quads = sizeof(LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim>>) * m_pLUT_lookup->get_num_quads();
                double byte_per_property = sizeof(double)*m_pLUT_lookup->m_props_unique_points_leaves.num_points;
                double byte_total = (byte_forest_leaves + byte_forest_nonleaves + byte_quads + byte_per_property*m_pLUT_lookup->m_props_unique_points_leaves.num_props);
                headinfo.memory_total = m_pLUT_lookup->byte2string(byte_total);
                headinfo.memory_leaves = m_pLUT_lookup->byte2string(byte_forest_leaves);
                headinfo.memory_nonLeaves = m_pLUT_lookup->byte2string(byte_forest_nonleaves);
                headinfo.memory_quads = m_pLUT_lookup->byte2string(byte_quads);
                headinfo.memory_properties = m_pLUT_lookup->byte2string(byte_per_property);
            }
                break;
            case 3:
            {
                LOOKUPTABLE_FOREST::LookUpTableForest_3D* m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_3D(file_lut, NULL,false));
                if(printSummary)m_pLUT_lookup->print_summary();
            }
                break;
            default:
            ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(file_lut));
                break;
        }
        // construct head info


        return headinfo;
    }

    ThermodynamicPropertiesVector cxThermal::UpdateState_TPX(const std::vector<double>& T, const std::vector<double>& p, const std::vector<double>& X, bool isMeshGrid)
    {
        ThermodynamicPropertiesVector stateVector;
        stateVector.fluidName = name();
        if (!isMeshGrid)
        {
            if (T.size() != p.size() && T.size() != X.size()) ERROR("The size of input T,p,X vectors are not identical in cxThermal::UpdateState_TPX(const std::vector<double>& T, const std::vector<double>& p, const std::vector<double>& X, bool isMeshGrid), can not process vector calculation, please check.");
            size_t N = T.size();
            stateVector.resize(N);
            MultiProgressBar multiBar(N);
            ThermodynamicProperties props;
#ifdef USE_OMP
            if(get_num_threads()>1)STATUS("Parallel computing, threads number: "<<get_num_threads());
#pragma omp parallel for private(props) shared(stateVector, T, p, X)
#endif
            for (size_t i = 0; i < N; ++i) {
                stateVector.T[i] = T[i];    stateVector.p[i] = p[i];   stateVector.X[i] = X[i];
                UpdateState_TPX(props, T[i], p[i], X[i]);
                stateVector.fill(props, i);
                if(m_isShowProgressBar)
                {
#ifdef USE_OMP
#pragma omp critical
#endif
                    multiBar.Update();
                }
            }
        }
        else
        {
            size_t nT = T.size(), nP = p.size(), nX = X.size();
            size_t N = nT*nP*nX;
            size_t nXnT = nT*nX;
            stateVector.resize(N);
            ThermodynamicProperties props;
            // size_t  ind = 0;
            MultiProgressBar multiBar(nT*nP);
#ifdef USE_OMP
            if(get_num_threads()>1)STATUS("Parallel computing, threads number: "<<get_num_threads()<<"\n");
#pragma omp parallel for private(props) shared(stateVector, T, p, X)
#endif
            for (int i = 0; i < nP; ++i) {
                for (int j = 0; j < nT; ++j) {
                    for (int k = 0; k < nX; ++k) {
                        size_t  ind = k + j*nX + i*(nXnT);
                        stateVector.T[ind] = T[j];    stateVector.p[ind] = p[i];   stateVector.X[ind] = X[k];
                        UpdateState_TPX(props, T[j], p[i], X[k]);
                        stateVector.fill(props, ind);
                    }
                    if(m_isShowProgressBar)
                    {
#ifdef USE_OMP
#pragma omp critical
#endif
                        multiBar.Update();
                    }
                }
            }
        }
        return stateVector;
    }

    ThermodynamicPropertiesVector cxThermal::UpdateState_HPX(const std::vector<double>& H, const std::vector<double>& p, const std::vector<double>& X, bool isMeshGrid)
    {
        ThermodynamicPropertiesVector stateVector;
        stateVector.fluidName = name();
        if (!isMeshGrid)
        {
            size_t N = H.size();
            stateVector.resize(N);
            MultiProgressBar multiBar(N);
            ThermodynamicProperties props;
#ifdef USE_OMP
            if(get_num_threads()>1)STATUS("Parallel computing, threads number: "<<get_num_threads());
#pragma omp parallel for private(props) shared(stateVector, H, p, X)
#endif
            for (size_t i = 0; i < N; ++i) {
                stateVector.H[i] = H[i];    stateVector.p[i] = p[i];   stateVector.X[i] = X[i];
                // 1. calculate phase region
                UpdateState_HPX(props, H[i], p[i], X[i]);
                stateVector.fill(props, i);
                if(m_isShowProgressBar)
                {
#ifdef USE_OMP
#pragma omp critical
#endif
                    multiBar.Update();
                }
            }
        } else
        {
            size_t nH = H.size(), nP = p.size(), nX = X.size();
            size_t N = nH*nP*nX;
            size_t nXnH = nX*nH;
            stateVector.resize(N);
            ThermodynamicProperties props;
            // size_t  ind = 0;
            MultiProgressBar multiBar(nH*nP);
#ifdef USE_OMP
            if(get_num_threads()>1)STATUS("Parallel computing, threads number: "<<get_num_threads()<<"\n");
#pragma omp parallel for private(props) shared(stateVector, H, p, X)
#endif
            for (int i = 0; i < nP; ++i) {
                for (int j = 0; j < nH; ++j) {
                    for (int k = 0; k < nX; ++k) {
                        size_t  ind = k + j*nX + i*(nXnH);
                        stateVector.H[ind] = H[j];    stateVector.p[ind] = p[i];   stateVector.X[ind] = X[k];
                        UpdateState_HPX(props, H[j], p[i], X[k]);
                        stateVector.fill(props, ind);
                    }
                    if(m_isShowProgressBar)
                    {
#ifdef USE_OMP
#pragma omp critical
#endif
                        multiBar.Update();
                    }
                }
            }
        }
        return stateVector;
    }

    void cxThermal::createLUT_3D(double *xyz_min, double *xyz_max, LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level,
                                int max_level, int update_which_props)
    {
        // parsing which properties need to be calculated
        parse_update_which_props(update_which_props);

        destroyLUT(m_pLUT, m_dim_lut); //destroy LUT pointer and release all related data if it exists.
        clock_t start = clock();
        STATUS("Creating 3D lookup table ...");
        m_dim_lut = 3;
        LOOKUPTABLE_FOREST::LookUpTableForest_3D* tmp_lut_3D = new LOOKUPTABLE_FOREST::LookUpTableForest_3D (xyz_min, xyz_max, TorH, max_level, m_update_which_props, this);
        m_pLUT = tmp_lut_3D;
        // refine
        tmp_lut_3D->set_min_level(min_level);
        tmp_lut_3D->refine(refine_uniform);
        // parallel refine
        if(tmp_lut_3D->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
        {
#ifdef USE_OMP
#pragma omp parallel //shared(n)
#endif
            {
#ifdef USE_OMP
#pragma omp single
#endif
                {
                    printf("Do refinement using %d threads.\n", m_num_threads);
                    tmp_lut_3D->refine(RefineFunc_PTX);
                }
            }
            //
            STATUS_time("Lookup table refinement done", (clock() - start)/m_num_threads);
            tmp_lut_3D->construct_props_leaves(cal_prop_PTX);
        }else if (tmp_lut_3D->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_H)
        {
#ifdef USE_OMP
#pragma omp parallel //shared(n)
#endif
            {
#ifdef USE_OMP
#pragma omp single
#endif
                {
                    printf("Do refinement using %d threads.\n", m_num_threads);
                    tmp_lut_3D->refine(RefineFunc_PHX);
                }
            }
            //
            STATUS_time("Lookup table refinement done", (clock() - start)/m_num_threads);
            tmp_lut_3D->construct_props_leaves(cal_prop_PHX);
        }else
        {
            ERROR("The EOS space only support TPX and HPX!");
        }
        tmp_lut_3D->print_summary();
    }

    void cxThermal::createLUT_3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                                LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level, int max_level,
                                int update_which_props)
    {
        double xyz_min[3] = {xmin, ymin, zmin};
        double xyz_max[3] = {xmax, ymax, zmax};
        createLUT_3D(xyz_min, xyz_max, TorH,  min_level, max_level, update_which_props);
    }

    void cxThermal::writeMeshGrid2VTK(const string &vtkFile,
                                     const vector<double> &x, const std::string& xTitle,
                                     const vector<double> &y, const std::string& yTitle,
                                     const vector<double> &z, const std::string& zTitle,
                                     const vector<std::vector<double>> &props, const vector<propInfo> &propsInfo,
                                     bool isNormalize) {
        using namespace std;
        STATUS("Writing mesh grid to structured vtk grid file : "<<vtkFile);
        if (props.empty()) ERROR("There is no properties for mesh grid writing!");
        ofstream fpout(vtkFile);
        if(!fpout)
        {
            ERROR("Can not open file: "<<vtkFile);
            exit(0);
        }
        string fname_py=vtkFile+".py";
        //  write vtk head
        fpout<<"# vtk DataFile Version 2.0"<<endl;
        fpout<<"Properties of seawater"<<endl;
        fpout<<"ASCII"<<endl;
        fpout<<"DATASET RECTILINEAR_GRID"<<endl;
        fpout<<"DIMENSIONS "<<x.size()<<" "<<y.size()<<" "<<z.size()<<endl;
        double len_x=1, len_y=1, len_z=1;
        double xMAX=*std::max_element(x.begin(), x.end());
        double xMIN=*std::min_element(x.begin(), x.end());
        double yMAX=*std::max_element(y.begin(), y.end());
        double yMIN=*std::min_element(y.begin(), y.end());
        double zMAX=*std::max_element(z.begin(), z.end());
        double zMIN=*std::min_element(z.begin(), z.end());
        double scale_x=1, scale_y=1, scale_z=1;
        len_x=(xMAX==xMIN ? 1: xMAX-xMIN);
        len_y=(yMAX==yMIN ? 1: yMAX-yMIN);
        len_z=(zMAX==zMIN ? 1: zMAX-zMIN);
        if(!isNormalize) // if not normalize the vtk file, write a python script to better visualize result
        {
            scale_y=len_x/len_y;
            scale_z=len_x/len_z;
            ofstream fout_py(fname_py);
            if(!fout_py)
            {
                cout<<"Warning: cannot generate pvPython script for Paraview. "<<fname_py<<endl;
            }else
            {
                fout_py<<"from paraview.simple import *"<<endl;
                fout_py<<"xHvtk = LegacyVTKReader(FileNames=[\'"<<vtkFile<<"\'])"<<endl;
                fout_py<<"renderView1 = GetActiveViewOrCreate(\'RenderView\')"<<endl;
                fout_py<<"xHvtkDisplay = Show(xHvtk, renderView1)"<<endl;
                fout_py<<"xHvtkDisplay.Representation = \'Surface\'"<<endl;
                fout_py<<"renderView1.AxesGrid.Visibility = 1"<<endl;
                fout_py<<"xHvtkDisplay.Scale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                fout_py<<"renderView1.AxesGrid.DataScale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                // fout_py<<"renderView1.AxesGrid.DataBoundsInflateFactor = 0"<<endl;
                fout_py<<"renderView1.AxesGrid.XTitle = \'"<<xTitle<<"\'"<<endl;
                fout_py<<"renderView1.AxesGrid.YTitle = \'"<<yTitle<<"\'"<<endl;
                fout_py<<"renderView1.AxesGrid.ZTitle = \'"<<zTitle<<"\'"<<endl;
                if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleFontSize = 16"<<endl;
                if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleBold = 1"<<endl;
                if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleFontSize = 16"<<endl;
                if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleBold = 1"<<endl;
                if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleFontSize = 16"<<endl;
                if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleBold = 1"<<endl;
                // set default data source as the first property
                fout_py<<"#set default data source as "<<propsInfo[0].shortName<<propsInfo[0].unit<<endl;
                fout_py<<"paraview.simple._DisableFirstRenderCameraReset()"<<endl;
                fout_py<<"legacyVTKReader1 = GetActiveSource()"<<endl;
                fout_py<<"renderView1 = GetActiveViewOrCreate('RenderView')"<<endl;
                fout_py<<"legacyVTKReader1Display = GetDisplayProperties(legacyVTKReader1, view=renderView1)"<<endl;
                fout_py<<"ColorBy(legacyVTKReader1Display, ('POINTS', '"<<propsInfo[0].shortName<<propsInfo[0].unit<<"'))"<<endl;
                fout_py<<"legacyVTKReader1Display.RescaleTransferFunctionToDataRange(True, False)"<<endl;
                fout_py<<"legacyVTKReader1Display.SetScalarBarVisibility(renderView1, True)"<<endl;
                // fout_py<<"phaseRegionLUT = GetColorTransferFunction('"<<propsInfo[0].shortName<<propsInfo[0].unit<<"')\n"<<endl;
                fout_py<<"renderView1.ResetCamera()"<<endl;
                fout_py.close();
                STATUS("Paraview-python script is generated as : "<<fname_py);
            }
        }
        fpout<<"X_COORDINATES "<<x.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(double i : x)fpout<<(i-xMIN)/len_x<<" ";fpout<<endl;
        }else
        {
            for(double i : x)fpout<<i<<" ";fpout<<endl;
        }
        fpout<<"Y_COORDINATES "<<y.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(double i : y)fpout<<(i-yMIN)/len_y<<" ";fpout<<endl;
        }else
        {
            for(double i : y)fpout<<i<<" ";fpout<<endl;
        }
        fpout<<"Z_COORDINATES "<<z.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(double i : z)fpout<<(i-zMIN)/len_z<<" ";fpout<<endl;
        }else
        {
            for(double i : z)fpout<<i<<" ";fpout<<endl;
        }
        fpout<<"POINT_DATA "<<props.size()<<endl;
        // 1.for loop: write point data of properties
        for (int i = 0; i < propsInfo.size(); ++i) {
            fpout<<"SCALARS "<<propsInfo[i].shortName<<propsInfo[i].unit<<" double"<<endl;
            fpout<<"LOOKUP_TABLE default"<<endl;
            for (const auto & prop : props) {
                fpout<<prop[i]<<" ";
            }
            fpout<<endl;
        }
        fpout.close();
        if(!isNormalize)STATUS("You can use command of "<<COLOR_PURPLE<<"paraview --script="<<fname_py<<COLOR_DEFAULT<<" to visualize result in paraview")

    }

    void cxThermal::UpdateState_TPX(ThermodynamicPropertiesArray &stateArray, const size_t &N, const double *T, const double *p, const double *X)
    {
        stateArray.fluidName = name();
        stateArray.N = N;
        ThermodynamicProperties props;
        for (int i = 0; i < N; ++i) {
            UpdateState_TPX(props, T[i], p[i], X[i]);
            stateArray.fill(props, i);
        }

    }

    void cxThermal::UpdateState_HPX(ThermodynamicPropertiesArray &stateArray, const size_t &N, const double *H, const double *p, const double *X)
    {
        stateArray.fluidName = name();
        stateArray.N = N;
        ThermodynamicProperties props;
        for (int i = 0; i < N; ++i) {
            UpdateState_HPX(props, H[i], p[i], X[i]);
            stateArray.fill(props, i);
        }
    }

    INDEX_FLUID cxThermal::validateFluid(const std::string& fluidName) {
        std::map<std::string, INDEX_FLUID> supportedFluids;
        // ------- initialize supported fluids ---------
        supportedFluids["Water"] = FLUID_Water;
        supportedFluids["H2O-NaCl"] = FLUID_H2O_NaCl;
        // for other fluids in CoolProp, see http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
        //----------------------------------------------
        if(supportedFluids.count(fluidName)!=0)
        {
            return supportedFluids[fluidName];
        }
        return Fluid_Unknown;
    }

    void *cxThermal::get_pLUT() {
        if (!m_pLUT)
        {
            WARNING("The member variable m_pLUT of thermo class is NULL, are you sure the calling function of get_pLUT is a proper one?\nIf you want to get pointer of LUT for property lookup, please use get_pLUT_lookup()!")
        }
        return m_pLUT;
    }

    void *cxThermal::get_pLUT_lookup() {
        if(!m_pLUT_lookup)
        {
            WARNING("The member variable m_pLUT_lookup of thermo class is NULL, are you sure the calling function of get_pLUT_lookup is a proper one?\nThis function is usually used for properties look up from a LUT binary file, do you mean another function get_pLUT()?")
        }
        return m_pLUT_lookup;
    }

}