

# File triangle.cpp



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**triangle.cpp**](triangle_8cpp.md)

[Go to the source code of this file](triangle_8cpp_source.md)



* `#include <stdio.h>`
* `#include <stdlib.h>`
* `#include <string.h>`
* `#include <math.h>`
* `#include <sys/time.h>`
* `#include "triangle.h"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**badsubseg**](structbadsubseg.md) <br> |
| struct | [**badtriang**](structbadtriang.md) <br> |
| struct | [**behavior**](structbehavior.md) <br> |
| struct | [**event**](structevent.md) <br> |
| struct | [**flipstacker**](structflipstacker.md) <br> |
| struct | [**memorypool**](structmemorypool.md) <br> |
| struct | [**mesh**](structmesh.md) <br> |
| struct | [**osub**](structosub.md) <br> |
| struct | [**otri**](structotri.md) <br> |
| struct | [**splaynode**](structsplaynode.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**finddirectionresult**](#enum-finddirectionresult)  <br> |
| enum  | [**insertvertexresult**](#enum-insertvertexresult)  <br> |
| enum  | [**locateresult**](#enum-locateresult)  <br> |
| typedef REAL \*\* | [**subseg**](#typedef-subseg)  <br> |
| typedef REAL \*\* | [**triangle**](#typedef-triangle)  <br> |
| typedef REAL \* | [**vertex**](#typedef-vertex)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  REAL | [**ccwerrboundA**](#variable-ccwerrbounda)  <br> |
|  REAL | [**ccwerrboundB**](#variable-ccwerrboundb)  <br> |
|  REAL | [**ccwerrboundC**](#variable-ccwerrboundc)  <br> |
|  REAL | [**epsilon**](#variable-epsilon)  <br> |
|  REAL | [**iccerrboundA**](#variable-iccerrbounda)  <br> |
|  REAL | [**iccerrboundB**](#variable-iccerrboundb)  <br> |
|  REAL | [**iccerrboundC**](#variable-iccerrboundc)  <br> |
|  int | [**minus1mod3**](#variable-minus1mod3)   = `{2, 0, 1}`<br> |
|  REAL | [**o3derrboundA**](#variable-o3derrbounda)  <br> |
|  REAL | [**o3derrboundB**](#variable-o3derrboundb)  <br> |
|  REAL | [**o3derrboundC**](#variable-o3derrboundc)  <br> |
|  int | [**plus1mod3**](#variable-plus1mod3)   = `{1, 2, 0}`<br> |
|  unsigned long | [**randomseed**](#variable-randomseed)  <br> |
|  REAL | [**resulterrbound**](#variable-resulterrbound)  <br> |
|  REAL | [**splitter**](#variable-splitter)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**alternateaxes**](#function-alternateaxes) (vertex \* sortarray, int arraysize, int axis) <br> |
|  void | [**badsubsegdealloc**](#function-badsubsegdealloc) (struct [**mesh**](structmesh.md) \* m, struct [**badsubseg**](structbadsubseg.md) \* dyingseg) <br> |
|  struct [**badsubseg**](structbadsubseg.md) \* | [**badsubsegtraverse**](#function-badsubsegtraverse) (struct [**mesh**](structmesh.md) \* m) <br> |
|  void | [**boundingbox**](#function-boundingbox) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**carveholes**](#function-carveholes) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, REAL \* holelist, int holes, REAL \* regionlist, int regions) <br> |
|  void | [**check4deadevent**](#function-check4deadevent) (struct [**otri**](structotri.md) \* checktri, struct [**event**](structevent.md) \*\* freeevents, struct [**event**](structevent.md) \*\* eventheap, int \* heapsize) <br> |
|  void | [**checkdelaunay**](#function-checkdelaunay) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**checkmesh**](#function-checkmesh) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  int | [**checkseg4encroach**](#function-checkseg4encroach) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**osub**](structosub.md) \* testsubseg) <br> |
|  REAL | [**circletop**](#function-circletop) (struct [**mesh**](structmesh.md) \* m, vertex pa, vertex pb, vertex pc, REAL ccwabc) <br> |
|  struct [**splaynode**](structsplaynode.md) \* | [**circletopinsert**](#function-circletopinsert) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**splaynode**](structsplaynode.md) \* splayroot, struct [**otri**](structotri.md) \* newkey, vertex pa, vertex pb, vertex pc, REAL topy) <br> |
|  void | [**conformingedge**](#function-conformingedge) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex endpoint1, vertex endpoint2, int newmark) <br> |
|  void | [**constrainededge**](#function-constrainededge) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* starttri, vertex endpoint2, int newmark) <br> |
|  REAL | [**counterclockwise**](#function-counterclockwise) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex pa, vertex pb, vertex pc) <br> |
|  REAL | [**counterclockwiseadapt**](#function-counterclockwiseadapt) (vertex pa, vertex pb, vertex pc, REAL detsum) <br> |
|  void | [**createeventheap**](#function-createeventheap) (struct [**mesh**](structmesh.md) \* m, struct [**event**](structevent.md) \*\*\* eventheap, struct [**event**](structevent.md) \*\* events, struct [**event**](structevent.md) \*\* freeevents) <br> |
|  long | [**delaunay**](#function-delaunay) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**delaunayfixup**](#function-delaunayfixup) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* fixuptri, int leftside) <br> |
|  void | [**deletevertex**](#function-deletevertex) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* deltri) <br> |
|  struct [**badtriang**](structbadtriang.md) \* | [**dequeuebadtriang**](#function-dequeuebadtriang) (struct [**mesh**](structmesh.md) \* m) <br> |
|  long | [**divconqdelaunay**](#function-divconqdelaunay) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**divconqrecurse**](#function-divconqrecurse) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex \* sortarray, int vertices, int axis, struct [**otri**](structotri.md) \* farleft, struct [**otri**](structotri.md) \* farright) <br> |
|  void | [**dummyinit**](#function-dummyinit) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int trianglebytes, int subsegbytes) <br> |
|  void | [**enforcequality**](#function-enforcequality) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**enqueuebadtri**](#function-enqueuebadtri) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* enqtri, REAL minedge, vertex enqapex, vertex enqorg, vertex enqdest) <br> |
|  void | [**enqueuebadtriang**](#function-enqueuebadtriang) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**badtriang**](structbadtriang.md) \* badtri) <br> |
|  REAL | [**estimate**](#function-estimate) (int elen, REAL \* e) <br> |
|  void | [**eventheapdelete**](#function-eventheapdelete) (struct [**event**](structevent.md) \*\* heap, int heapsize, int eventnum) <br> |
|  void | [**eventheapify**](#function-eventheapify) (struct [**event**](structevent.md) \*\* heap, int heapsize, int eventnum) <br> |
|  void | [**eventheapinsert**](#function-eventheapinsert) (struct [**event**](structevent.md) \*\* heap, int heapsize, struct [**event**](structevent.md) \* newevent) <br> |
|  void | [**exactinit**](#function-exactinit) () <br> |
|  int | [**fast\_expansion\_sum\_zeroelim**](#function-fast_expansion_sum_zeroelim) (int elen, REAL \* e, int flen, REAL \* f, REAL \* h) <br> |
|  void | [**findcircumcenter**](#function-findcircumcenter) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex torg, vertex tdest, vertex tapex, vertex circumcenter, REAL \* xi, REAL \* eta, int offcenter) <br> |
|  enum finddirectionresult | [**finddirection**](#function-finddirection) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* searchtri, vertex searchpoint) <br> |
|  void | [**flip**](#function-flip) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* flipedge) <br> |
|  void | [**formskeleton**](#function-formskeleton) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \* segmentlist, int \* segmentmarkerlist, int numberofsegments) <br> |
|  struct [**splaynode**](structsplaynode.md) \* | [**frontlocate**](#function-frontlocate) (struct [**mesh**](structmesh.md) \* m, struct [**splaynode**](structsplaynode.md) \* splayroot, struct [**otri**](structotri.md) \* bottommost, vertex searchvertex, struct [**otri**](structotri.md) \* searchtri, int \* farright) <br> |
|  vertex | [**getvertex**](#function-getvertex) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int number) <br> |
|  void | [**highorder**](#function-highorder) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  REAL | [**incircle**](#function-incircle) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex pa, vertex pb, vertex pc, vertex pd) <br> |
|  REAL | [**incircleadapt**](#function-incircleadapt) (vertex pa, vertex pb, vertex pc, vertex pd, REAL permanent) <br> |
|  long | [**incrementaldelaunay**](#function-incrementaldelaunay) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**infecthull**](#function-infecthull) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**initializetrisubpools**](#function-initializetrisubpools) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**initializevertexpool**](#function-initializevertexpool) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**insertsegment**](#function-insertsegment) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex endpoint1, vertex endpoint2, int newmark) <br> |
|  void | [**insertsubseg**](#function-insertsubseg) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* tri, int subsegmark) <br> |
|  enum insertvertexresult | [**insertvertex**](#function-insertvertex) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex newvertex, struct [**otri**](structotri.md) \* searchtri, struct [**osub**](structosub.md) \* splitseg, int segmentflaws, int triflaws) <br> |
|  void | [**internalerror**](#function-internalerror) () <br> |
|  enum locateresult | [**locate**](#function-locate) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex searchpoint, struct [**otri**](structotri.md) \* searchtri) <br> |
|  void | [**makesubseg**](#function-makesubseg) (struct [**mesh**](structmesh.md) \* m, struct [**osub**](structosub.md) \* newsubseg) <br> |
|  void | [**maketriangle**](#function-maketriangle) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* newotri) <br> |
|  void | [**makevertexmap**](#function-makevertexmap) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**markhull**](#function-markhull) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**mergehulls**](#function-mergehulls) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* farleft, struct [**otri**](structotri.md) \* innerleft, struct [**otri**](structotri.md) \* innerright, struct [**otri**](structotri.md) \* farright, int axis) <br> |
|  REAL | [**nonregular**](#function-nonregular) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex pa, vertex pb, vertex pc, vertex pd) <br> |
|  void | [**numbernodes**](#function-numbernodes) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  REAL | [**orient3d**](#function-orient3d) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex pa, vertex pb, vertex pc, vertex pd, REAL aheight, REAL bheight, REAL cheight, REAL dheight) <br> |
|  REAL | [**orient3dadapt**](#function-orient3dadapt) (vertex pa, vertex pb, vertex pc, vertex pd, REAL aheight, REAL bheight, REAL cheight, REAL dheight, REAL permanent) <br> |
|  void | [**parsecommandline**](#function-parsecommandline) (int argc, char \*\* argv, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**plague**](#function-plague) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  VOID \* | [**poolalloc**](#function-poolalloc) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  void | [**pooldealloc**](#function-pooldealloc) (struct [**memorypool**](structmemorypool.md) \* pool, VOID \* dyingitem) <br> |
|  void | [**pooldeinit**](#function-pooldeinit) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  void | [**poolinit**](#function-poolinit) (struct [**memorypool**](structmemorypool.md) \* pool, int bytecount, int itemcount, int firstitemcount, int alignment) <br> |
|  void | [**poolrestart**](#function-poolrestart) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  void | [**poolzero**](#function-poolzero) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  enum locateresult | [**preciselocate**](#function-preciselocate) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, vertex searchpoint, struct [**otri**](structotri.md) \* searchtri, int stopatsubsegment) <br> |
|  void | [**precisionerror**](#function-precisionerror) () <br> |
|  void | [**printsubseg**](#function-printsubseg) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**osub**](structosub.md) \* s) <br> |
|  void | [**printtriangle**](#function-printtriangle) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* t) <br> |
|  void | [**quality\_statistics**](#function-quality_statistics) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  unsigned long | [**randomnation**](#function-randomnation) (unsigned int choices) <br> |
|  int | [**reconstruct**](#function-reconstruct) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \* trianglelist, REAL \* triangleattriblist, REAL \* trianglearealist, int elements, int corners, int attribs, int \* segmentlist, int \* segmentmarkerlist, int numberofsegments) <br> |
|  void | [**regionplague**](#function-regionplague) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, REAL attribute, REAL area) <br> |
|  long | [**removebox**](#function-removebox) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  long | [**removeghosts**](#function-removeghosts) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* startghost) <br> |
|  void | [**report**](#function-report) (struct [**triangulateio**](structtriangulateio.md) \* io, int markers, int reporttriangles, int reportneighbors, int reportsegments, int reportedges, int reportnorms) <br> |
|  int | [**rightofhyperbola**](#function-rightofhyperbola) (struct [**mesh**](structmesh.md) \* m, struct [**otri**](structotri.md) \* fronttri, vertex newsite) <br> |
|  int | [**scale\_expansion\_zeroelim**](#function-scale_expansion_zeroelim) (int elen, REAL \* e, REAL b, REAL \* h) <br> |
|  int | [**scoutsegment**](#function-scoutsegment) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* searchtri, vertex endpoint2, int newmark) <br> |
|  void | [**segmentintersection**](#function-segmentintersection) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* splittri, struct [**osub**](structosub.md) \* splitsubseg, vertex endpoint2) <br> |
|  struct [**splaynode**](structsplaynode.md) \* | [**splay**](#function-splay) (struct [**mesh**](structmesh.md) \* m, struct [**splaynode**](structsplaynode.md) \* splaytree, vertex searchpoint, struct [**otri**](structotri.md) \* searchtri) <br> |
|  struct [**splaynode**](structsplaynode.md) \* | [**splayinsert**](#function-splayinsert) (struct [**mesh**](structmesh.md) \* m, struct [**splaynode**](structsplaynode.md) \* splayroot, struct [**otri**](structotri.md) \* newkey, vertex searchpoint) <br> |
|  void | [**splitencsegs**](#function-splitencsegs) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int triflaws) <br> |
|  void | [**splittriangle**](#function-splittriangle) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**badtriang**](structbadtriang.md) \* badtri) <br> |
|  void | [**statistics**](#function-statistics) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**subsegdealloc**](#function-subsegdealloc) (struct [**mesh**](structmesh.md) \* m, subseg \* dyingsubseg) <br> |
|  subseg \* | [**subsegtraverse**](#function-subsegtraverse) (struct [**mesh**](structmesh.md) \* m) <br> |
|  long | [**sweeplinedelaunay**](#function-sweeplinedelaunay) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**tallyencs**](#function-tallyencs) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**tallyfaces**](#function-tallyfaces) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**testtriangle**](#function-testtriangle) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* testtri) <br> |
|  void | [**transfernodes**](#function-transfernodes) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, REAL \* pointlist, REAL \* pointattriblist, int \* pointmarkerlist, int numberofpoints, int numberofpointattribs) <br> |
|  void | [**traversalinit**](#function-traversalinit) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  VOID \* | [**traverse**](#function-traverse) (struct [**memorypool**](structmemorypool.md) \* pool) <br> |
|  void | [**triangledealloc**](#function-triangledealloc) (struct [**mesh**](structmesh.md) \* m, triangle \* dyingtriangle) <br> |
|  void | [**triangledeinit**](#function-triangledeinit) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**triangleinit**](#function-triangleinit) (struct [**mesh**](structmesh.md) \* m) <br> |
|  triangle \* | [**triangletraverse**](#function-triangletraverse) (struct [**mesh**](structmesh.md) \* m) <br> |
|  void | [**triangulate**](#function-triangulate) (char \* triswitches, struct [**triangulateio**](structtriangulateio.md) \* in, struct [**triangulateio**](structtriangulateio.md) \* out, struct [**triangulateio**](structtriangulateio.md) \* vorout) <br> |
|  void | [**triangulatepolygon**](#function-triangulatepolygon) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* firstedge, struct [**otri**](structotri.md) \* lastedge, int edgecount, int doflip, int triflaws) <br> |
|  void | [**triexit**](#function-triexit) (int status) <br> |
|  void | [**trifree**](#function-trifree) (VOID \* memptr) <br> |
|  VOID \* | [**trimalloc**](#function-trimalloc) (int size) <br> |
|  int | [**triunsuitable**](#function-triunsuitable) (vertex triorg, vertex tridest, vertex triapex, REAL area) <br> |
|  void | [**undovertex**](#function-undovertex) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b) <br> |
|  void | [**unflip**](#function-unflip) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, struct [**otri**](structotri.md) \* flipedge) <br> |
|  void | [**vertexdealloc**](#function-vertexdealloc) (struct [**mesh**](structmesh.md) \* m, vertex dyingvertex) <br> |
|  void | [**vertexmedian**](#function-vertexmedian) (vertex \* sortarray, int arraysize, int median, int axis) <br> |
|  void | [**vertexsort**](#function-vertexsort) (vertex \* sortarray, int arraysize) <br> |
|  vertex | [**vertextraverse**](#function-vertextraverse) (struct [**mesh**](structmesh.md) \* m) <br> |
|  void | [**writeedges**](#function-writeedges) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \*\* edgelist, int \*\* edgemarkerlist) <br> |
|  void | [**writeelements**](#function-writeelements) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \*\* trianglelist, REAL \*\* triangleattriblist) <br> |
|  void | [**writeneighbors**](#function-writeneighbors) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \*\* neighborlist) <br> |
|  void | [**writenodes**](#function-writenodes) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, REAL \*\* pointlist, REAL \*\* pointattriblist, int \*\* pointmarkerlist) <br> |
|  void | [**writepoly**](#function-writepoly) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, int \*\* segmentlist, int \*\* segmentmarkerlist) <br> |
|  void | [**writevoronoi**](#function-writevoronoi) (struct [**mesh**](structmesh.md) \* m, struct [**behavior**](structbehavior.md) \* b, REAL \*\* vpointlist, REAL \*\* vpointattriblist, int \*\* vpointmarkerlist, int \*\* vedgelist, int \*\* vedgemarkerlist, REAL \*\* vnormlist) <br> |



























## Macros

| Type | Name |
| ---: | :--- |
| define  | [**Absolute**](triangle_8cpp.md#define-absolute) (a) `((a) &gt;= 0.0 ? (a) : -(a))`<br> |
| define  | [**BADSUBSEGPERBLOCK**](triangle_8cpp.md#define-badsubsegperblock)  `252`<br> |
| define  | [**BADTRIPERBLOCK**](triangle_8cpp.md#define-badtriperblock)  `4092`<br> |
| define  | [**DEADVERTEX**](triangle_8cpp.md#define-deadvertex)  `-32768`<br> |
| define  | [**FILENAMESIZE**](triangle_8cpp.md#define-filenamesize)  `2048`<br> |
| define  | [**FLIPSTACKERPERBLOCK**](triangle_8cpp.md#define-flipstackerperblock)  `252`<br> |
| define  | [**FREEVERTEX**](triangle_8cpp.md#define-freevertex)  `2`<br> |
| define  | [**Fast\_Two\_Sum**](triangle_8cpp.md#define-fast_two_sum) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Fast\_Two\_Sum\_Tail**](triangle_8cpp.md#define-fast_two_sum_tail) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**INEXACT**](triangle_8cpp.md#define-inexact)  `/\* Nothing \*/`<br> |
| define  | [**INPUTLINESIZE**](triangle_8cpp.md#define-inputlinesize)  `1024`<br> |
| define  | [**INPUTVERTEX**](triangle_8cpp.md#define-inputvertex)  `0`<br> |
| define  | [**ONETHIRD**](triangle_8cpp.md#define-onethird)  `0.333333333333333333333333333333333333333333333333333333333333`<br> |
| define  | [**PI**](triangle_8cpp.md#define-pi)  `3.141592653589793238462643383279502884197169399375105820974944592308`<br> |
| define  | [**REAL**](triangle_8cpp.md#define-real)  `double`<br> |
| define  | [**SAMPLEFACTOR**](triangle_8cpp.md#define-samplefactor)  `11`<br> |
| define  | [**SAMPLERATE**](triangle_8cpp.md#define-samplerate)  `10`<br> |
| define  | [**SEGMENTVERTEX**](triangle_8cpp.md#define-segmentvertex)  `1`<br> |
| define  | [**SPLAYNODEPERBLOCK**](triangle_8cpp.md#define-splaynodeperblock)  `508`<br> |
| define  | [**SQUAREROOTTWO**](triangle_8cpp.md#define-squareroottwo)  `1.4142135623730950488016887242096980785696718753769480732`<br> |
| define  | [**STARTINDEX**](triangle_8cpp.md#define-startindex)  `0`<br> |
| define  | [**SUBSEGPERBLOCK**](triangle_8cpp.md#define-subsegperblock)  `508       /\* Number of subsegments allocated at once. \*/`<br> |
| define  | [**Split**](triangle_8cpp.md#define-split) (a, ahi, alo) `/* multi line expression */`<br> |
| define  | [**Square**](triangle_8cpp.md#define-square) (a, x, y) `/* multi line expression */`<br> |
| define  | [**Square\_Tail**](triangle_8cpp.md#define-square_tail) (a, x, y) `/* multi line expression */`<br> |
| define  | [**TRILIBRARY**](triangle_8cpp.md#define-trilibrary)  `1`<br> |
| define  | [**TRIPERBLOCK**](triangle_8cpp.md#define-triperblock)  `4092           /\* Number of triangles allocated at once. \*/`<br> |
| define  | [**Two\_Diff**](triangle_8cpp.md#define-two_diff) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Diff\_Tail**](triangle_8cpp.md#define-two_diff_tail) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_One\_Diff**](triangle_8cpp.md#define-two_one_diff) (a1, a0, b, x2, x1, x0) `/* multi line expression */`<br> |
| define  | [**Two\_One\_Product**](triangle_8cpp.md#define-two_one_product) (a1, a0, b, x3, x2, x1, x0) `/* multi line expression */`<br> |
| define  | [**Two\_One\_Sum**](triangle_8cpp.md#define-two_one_sum) (a1, a0, b, x2, x1, x0) `/* multi line expression */`<br> |
| define  | [**Two\_Product**](triangle_8cpp.md#define-two_product) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Product\_Presplit**](triangle_8cpp.md#define-two_product_presplit) (a, b, bhi, blo, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Product\_Tail**](triangle_8cpp.md#define-two_product_tail) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Sum**](triangle_8cpp.md#define-two_sum) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Sum\_Tail**](triangle_8cpp.md#define-two_sum_tail) (a, b, x, y) `/* multi line expression */`<br> |
| define  | [**Two\_Two\_Diff**](triangle_8cpp.md#define-two_two_diff) (a1, a0, b1, b0, x3, x2, x1, x0) `/* multi line expression */`<br> |
| define  | [**Two\_Two\_Sum**](triangle_8cpp.md#define-two_two_sum) (a1, a0, b1, b0, x3, x2, x1, x0) `/* multi line expression */`<br> |
| define  | [**UNDEADVERTEX**](triangle_8cpp.md#define-undeadvertex)  `-32767`<br> |
| define  | [**VERTEXPERBLOCK**](triangle_8cpp.md#define-vertexperblock)  `4092         /\* Number of vertices allocated at once. \*/`<br> |
| define  | [**VIRUSPERBLOCK**](triangle_8cpp.md#define-virusperblock)  `1020   /\* Number of virus triangles allocated at once. \*/`<br> |
| define  | [**VOID**](triangle_8cpp.md#define-void)  `int`<br> |
| define  | [**apex**](triangle_8cpp.md#define-apex) (otri, vertexptr) `vertexptr = (vertex) ([**otri**](structotri.md)).tri[([**otri**](structotri.md)).orient + 3]`<br> |
| define  | [**areabound**](triangle_8cpp.md#define-areabound) (otri) `((REAL \*) ([**otri**](structotri.md)).tri)[m-&gt;areaboundindex]`<br> |
| define  | [**bond**](triangle_8cpp.md#define-bond) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**deadsubseg**](triangle_8cpp.md#define-deadsubseg) (sub) `((sub)[1] == (subseg) NULL)`<br> |
| define  | [**deadtri**](triangle_8cpp.md#define-deadtri) (tria) `((tria)[1] == (triangle) NULL)`<br> |
| define  | [**decode**](triangle_8cpp.md#define-decode) (ptr, otri) `/* multi line expression */`<br> |
| define  | [**dest**](triangle_8cpp.md#define-dest) (otri, vertexptr) `vertexptr = (vertex) ([**otri**](structotri.md)).tri[minus1mod3[([**otri**](structotri.md)).orient] + 3]`<br> |
| define  | [**dissolve**](triangle_8cpp.md#define-dissolve) (otri) `([**otri**](structotri.md)).tri[([**otri**](structotri.md)).orient] = (triangle) m-&gt;dummytri`<br> |
| define  | [**dnext**](triangle_8cpp.md#define-dnext) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**dnextself**](triangle_8cpp.md#define-dnextself) (otri) `/* multi line expression */`<br> |
| define  | [**dprev**](triangle_8cpp.md#define-dprev) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**dprevself**](triangle_8cpp.md#define-dprevself) (otri) `/* multi line expression */`<br> |
| define  | [**elemattribute**](triangle_8cpp.md#define-elemattribute) (otri, attnum) `((REAL \*) ([**otri**](structotri.md)).tri)[m-&gt;elemattribindex + (attnum)]`<br> |
| define  | [**encode**](triangle_8cpp.md#define-encode) (otri) `(triangle) ((unsigned long) ([**otri**](structotri.md)).tri \| (unsigned long) ([**otri**](structotri.md)).orient)`<br> |
| define  | [**infect**](triangle_8cpp.md#define-infect) (otri) `/* multi line expression */`<br> |
| define  | [**infected**](triangle_8cpp.md#define-infected) (otri) `(((unsigned long) ([**otri**](structotri.md)).tri[6] & (unsigned long) 2l) != 0l)`<br> |
| define  | [**killsubseg**](triangle_8cpp.md#define-killsubseg) (sub) `/* multi line expression */`<br> |
| define  | [**killtri**](triangle_8cpp.md#define-killtri) (tria) `/* multi line expression */`<br> |
| define  | [**lnext**](triangle_8cpp.md#define-lnext) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**lnextself**](triangle_8cpp.md#define-lnextself) (otri) `([**otri**](structotri.md)).orient = plus1mod3[([**otri**](structotri.md)).orient]`<br> |
| define  | [**lprev**](triangle_8cpp.md#define-lprev) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**lprevself**](triangle_8cpp.md#define-lprevself) (otri) `([**otri**](structotri.md)).orient = minus1mod3[([**otri**](structotri.md)).orient]`<br> |
| define  | [**mark**](triangle_8cpp.md#define-mark) (osub) `(\* (int \*) (([**osub**](structosub.md)).ss + 8))`<br> |
| define  | [**onext**](triangle_8cpp.md#define-onext) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**onextself**](triangle_8cpp.md#define-onextself) (otri) `/* multi line expression */`<br> |
| define  | [**oprev**](triangle_8cpp.md#define-oprev) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**oprevself**](triangle_8cpp.md#define-oprevself) (otri) `/* multi line expression */`<br> |
| define  | [**org**](triangle_8cpp.md#define-org) (otri, vertexptr) `vertexptr = (vertex) ([**otri**](structotri.md)).tri[plus1mod3[([**otri**](structotri.md)).orient] + 3]`<br> |
| define  | [**otricopy**](triangle_8cpp.md#define-otricopy) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**otriequal**](triangle_8cpp.md#define-otriequal) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**rnext**](triangle_8cpp.md#define-rnext) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**rnextself**](triangle_8cpp.md#define-rnextself) (otri) `/* multi line expression */`<br> |
| define  | [**rprev**](triangle_8cpp.md#define-rprev) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**rprevself**](triangle_8cpp.md#define-rprevself) (otri) `/* multi line expression */`<br> |
| define  | [**sbond**](triangle_8cpp.md#define-sbond) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**sdecode**](triangle_8cpp.md#define-sdecode) (sptr, osub) `/* multi line expression */`<br> |
| define  | [**sdest**](triangle_8cpp.md#define-sdest) (osub, vertexptr) `vertexptr = (vertex) ([**osub**](structosub.md)).ss[3 - ([**osub**](structosub.md)).ssorient]`<br> |
| define  | [**sdissolve**](triangle_8cpp.md#define-sdissolve) (osub) `([**osub**](structosub.md)).ss[([**osub**](structosub.md)).ssorient] = (subseg) m-&gt;dummysub`<br> |
| define  | [**segdest**](triangle_8cpp.md#define-segdest) (osub, vertexptr) `vertexptr = (vertex) ([**osub**](structosub.md)).ss[5 - ([**osub**](structosub.md)).ssorient]`<br> |
| define  | [**segorg**](triangle_8cpp.md#define-segorg) (osub, vertexptr) `vertexptr = (vertex) ([**osub**](structosub.md)).ss[4 + ([**osub**](structosub.md)).ssorient]`<br> |
| define  | [**sencode**](triangle_8cpp.md#define-sencode) (osub) `(subseg) ((unsigned long) ([**osub**](structosub.md)).ss \| (unsigned long) ([**osub**](structosub.md)).ssorient)`<br> |
| define  | [**setapex**](triangle_8cpp.md#define-setapex) (otri, vertexptr) `([**otri**](structotri.md)).tri[([**otri**](structotri.md)).orient + 3] = (triangle) vertexptr`<br> |
| define  | [**setareabound**](triangle_8cpp.md#define-setareabound) (otri, value) `((REAL \*) ([**otri**](structotri.md)).tri)[m-&gt;areaboundindex] = value`<br> |
| define  | [**setdest**](triangle_8cpp.md#define-setdest) (otri, vertexptr) `([**otri**](structotri.md)).tri[minus1mod3[([**otri**](structotri.md)).orient] + 3] = (triangle) vertexptr`<br> |
| define  | [**setelemattribute**](triangle_8cpp.md#define-setelemattribute) (otri, attnum, value) `((REAL \*) ([**otri**](structotri.md)).tri)[m-&gt;elemattribindex + (attnum)] = value`<br> |
| define  | [**setmark**](triangle_8cpp.md#define-setmark) (osub, value) `\* (int \*) (([**osub**](structosub.md)).ss + 8) = value`<br> |
| define  | [**setorg**](triangle_8cpp.md#define-setorg) (otri, vertexptr) `([**otri**](structotri.md)).tri[plus1mod3[([**otri**](structotri.md)).orient] + 3] = (triangle) vertexptr`<br> |
| define  | [**setsdest**](triangle_8cpp.md#define-setsdest) (osub, vertexptr) `([**osub**](structosub.md)).ss[3 - ([**osub**](structosub.md)).ssorient] = (subseg) vertexptr`<br> |
| define  | [**setsegdest**](triangle_8cpp.md#define-setsegdest) (osub, vertexptr) `([**osub**](structosub.md)).ss[5 - ([**osub**](structosub.md)).ssorient] = (subseg) vertexptr`<br> |
| define  | [**setsegorg**](triangle_8cpp.md#define-setsegorg) (osub, vertexptr) `([**osub**](structosub.md)).ss[4 + ([**osub**](structosub.md)).ssorient] = (subseg) vertexptr`<br> |
| define  | [**setsorg**](triangle_8cpp.md#define-setsorg) (osub, vertexptr) `([**osub**](structosub.md)).ss[2 + ([**osub**](structosub.md)).ssorient] = (subseg) vertexptr`<br> |
| define  | [**setvertex2tri**](triangle_8cpp.md#define-setvertex2tri) (vx, value) `((triangle \*) (vx))[m-&gt;vertex2triindex] = value`<br> |
| define  | [**setvertexmark**](triangle_8cpp.md#define-setvertexmark) (vx, value) `((int \*) (vx))[m-&gt;vertexmarkindex] = value`<br> |
| define  | [**setvertextype**](triangle_8cpp.md#define-setvertextype) (vx, value) `((int \*) (vx))[m-&gt;vertexmarkindex + 1] = value`<br> |
| define  | [**snext**](triangle_8cpp.md#define-snext) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**snextself**](triangle_8cpp.md#define-snextself) (osub) `/* multi line expression */`<br> |
| define  | [**sorg**](triangle_8cpp.md#define-sorg) (osub, vertexptr) `vertexptr = (vertex) ([**osub**](structosub.md)).ss[2 + ([**osub**](structosub.md)).ssorient]`<br> |
| define  | [**spivot**](triangle_8cpp.md#define-spivot) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**spivotself**](triangle_8cpp.md#define-spivotself) (osub) `/* multi line expression */`<br> |
| define  | [**ssym**](triangle_8cpp.md#define-ssym) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**ssymself**](triangle_8cpp.md#define-ssymself) (osub) `([**osub**](structosub.md)).ssorient = 1 - ([**osub**](structosub.md)).ssorient`<br> |
| define  | [**stdissolve**](triangle_8cpp.md#define-stdissolve) (osub) `([**osub**](structosub.md)).ss[6 + ([**osub**](structosub.md)).ssorient] = (subseg) m-&gt;dummytri`<br> |
| define  | [**stpivot**](triangle_8cpp.md#define-stpivot) (osub, otri) `/* multi line expression */`<br> |
| define  | [**subsegcopy**](triangle_8cpp.md#define-subsegcopy) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**subsegequal**](triangle_8cpp.md#define-subsegequal) (osub1, osub2) `/* multi line expression */`<br> |
| define  | [**sym**](triangle_8cpp.md#define-sym) (otri1, otri2) `/* multi line expression */`<br> |
| define  | [**symself**](triangle_8cpp.md#define-symself) (otri) `/* multi line expression */`<br> |
| define  | [**tsbond**](triangle_8cpp.md#define-tsbond) (otri, osub) `/* multi line expression */`<br> |
| define  | [**tsdissolve**](triangle_8cpp.md#define-tsdissolve) (otri) `([**otri**](structotri.md)).tri[6 + ([**otri**](structotri.md)).orient] = (triangle) m-&gt;dummysub`<br> |
| define  | [**tspivot**](triangle_8cpp.md#define-tspivot) (otri, osub) `/* multi line expression */`<br> |
| define  | [**uninfect**](triangle_8cpp.md#define-uninfect) (otri) `/* multi line expression */`<br> |
| define  | [**vertex2tri**](triangle_8cpp.md#define-vertex2tri) (vx) `((triangle \*) (vx))[m-&gt;vertex2triindex]`<br> |
| define  | [**vertexmark**](triangle_8cpp.md#define-vertexmark) (vx) `((int \*) (vx))[m-&gt;vertexmarkindex]`<br> |
| define  | [**vertextype**](triangle_8cpp.md#define-vertextype) (vx) `((int \*) (vx))[m-&gt;vertexmarkindex + 1]`<br> |

## Public Types Documentation




### enum finddirectionresult 

```C++
enum finddirectionresult {
    WITHIN,
    LEFTCOLLINEAR,
    RIGHTCOLLINEAR
};
```




<hr>



### enum insertvertexresult 

```C++
enum insertvertexresult {
    SUCCESSFULVERTEX,
    ENCROACHINGVERTEX,
    VIOLATINGVERTEX,
    DUPLICATEVERTEX
};
```




<hr>



### enum locateresult 

```C++
enum locateresult {
    INTRIANGLE,
    ONEDGE,
    ONVERTEX,
    OUTSIDE
};
```




<hr>



### typedef subseg 

```C++
typedef REAL** subseg;
```




<hr>



### typedef triangle 

```C++
typedef REAL** triangle;
```




<hr>



### typedef vertex 

```C++
typedef REAL* vertex;
```




<hr>
## Public Attributes Documentation




### variable ccwerrboundA 

```C++
REAL ccwerrboundA;
```




<hr>



### variable ccwerrboundB 

```C++
REAL ccwerrboundB;
```




<hr>



### variable ccwerrboundC 

```C++
REAL ccwerrboundC;
```




<hr>



### variable epsilon 

```C++
REAL epsilon;
```




<hr>



### variable iccerrboundA 

```C++
REAL iccerrboundA;
```




<hr>



### variable iccerrboundB 

```C++
REAL iccerrboundB;
```




<hr>



### variable iccerrboundC 

```C++
REAL iccerrboundC;
```




<hr>



### variable minus1mod3 

```C++
int minus1mod3[3];
```




<hr>



### variable o3derrboundA 

```C++
REAL o3derrboundA;
```




<hr>



### variable o3derrboundB 

```C++
REAL o3derrboundB;
```




<hr>



### variable o3derrboundC 

```C++
REAL o3derrboundC;
```




<hr>



### variable plus1mod3 

```C++
int plus1mod3[3];
```




<hr>



### variable randomseed 

```C++
unsigned long randomseed;
```




<hr>



### variable resulterrbound 

```C++
REAL resulterrbound;
```




<hr>



### variable splitter 

```C++
REAL splitter;
```




<hr>
## Public Functions Documentation




### function alternateaxes 

```C++
void alternateaxes (
    vertex * sortarray,
    int arraysize,
    int axis
) 
```




<hr>



### function badsubsegdealloc 

```C++
void badsubsegdealloc (
    struct mesh * m,
    struct badsubseg * dyingseg
) 
```




<hr>



### function badsubsegtraverse 

```C++
struct badsubseg * badsubsegtraverse (
    struct mesh * m
) 
```




<hr>



### function boundingbox 

```C++
void boundingbox (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function carveholes 

```C++
void carveholes (
    struct mesh * m,
    struct behavior * b,
    REAL * holelist,
    int holes,
    REAL * regionlist,
    int regions
) 
```




<hr>



### function check4deadevent 

```C++
void check4deadevent (
    struct otri * checktri,
    struct event ** freeevents,
    struct event ** eventheap,
    int * heapsize
) 
```




<hr>



### function checkdelaunay 

```C++
void checkdelaunay (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function checkmesh 

```C++
void checkmesh (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function checkseg4encroach 

```C++
int checkseg4encroach (
    struct mesh * m,
    struct behavior * b,
    struct osub * testsubseg
) 
```




<hr>



### function circletop 

```C++
REAL circletop (
    struct mesh * m,
    vertex pa,
    vertex pb,
    vertex pc,
    REAL ccwabc
) 
```




<hr>



### function circletopinsert 

```C++
struct splaynode * circletopinsert (
    struct mesh * m,
    struct behavior * b,
    struct splaynode * splayroot,
    struct otri * newkey,
    vertex pa,
    vertex pb,
    vertex pc,
    REAL topy
) 
```




<hr>



### function conformingedge 

```C++
void conformingedge (
    struct mesh * m,
    struct behavior * b,
    vertex endpoint1,
    vertex endpoint2,
    int newmark
) 
```




<hr>



### function constrainededge 

```C++
void constrainededge (
    struct mesh * m,
    struct behavior * b,
    struct otri * starttri,
    vertex endpoint2,
    int newmark
) 
```




<hr>



### function counterclockwise 

```C++
REAL counterclockwise (
    struct mesh * m,
    struct behavior * b,
    vertex pa,
    vertex pb,
    vertex pc
) 
```




<hr>



### function counterclockwiseadapt 

```C++
REAL counterclockwiseadapt (
    vertex pa,
    vertex pb,
    vertex pc,
    REAL detsum
) 
```




<hr>



### function createeventheap 

```C++
void createeventheap (
    struct mesh * m,
    struct event *** eventheap,
    struct event ** events,
    struct event ** freeevents
) 
```




<hr>



### function delaunay 

```C++
long delaunay (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function delaunayfixup 

```C++
void delaunayfixup (
    struct mesh * m,
    struct behavior * b,
    struct otri * fixuptri,
    int leftside
) 
```




<hr>



### function deletevertex 

```C++
void deletevertex (
    struct mesh * m,
    struct behavior * b,
    struct otri * deltri
) 
```




<hr>



### function dequeuebadtriang 

```C++
struct badtriang * dequeuebadtriang (
    struct mesh * m
) 
```




<hr>



### function divconqdelaunay 

```C++
long divconqdelaunay (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function divconqrecurse 

```C++
void divconqrecurse (
    struct mesh * m,
    struct behavior * b,
    vertex * sortarray,
    int vertices,
    int axis,
    struct otri * farleft,
    struct otri * farright
) 
```




<hr>



### function dummyinit 

```C++
void dummyinit (
    struct mesh * m,
    struct behavior * b,
    int trianglebytes,
    int subsegbytes
) 
```




<hr>



### function enforcequality 

```C++
void enforcequality (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function enqueuebadtri 

```C++
void enqueuebadtri (
    struct mesh * m,
    struct behavior * b,
    struct otri * enqtri,
    REAL minedge,
    vertex enqapex,
    vertex enqorg,
    vertex enqdest
) 
```




<hr>



### function enqueuebadtriang 

```C++
void enqueuebadtriang (
    struct mesh * m,
    struct behavior * b,
    struct badtriang * badtri
) 
```




<hr>



### function estimate 

```C++
REAL estimate (
    int elen,
    REAL * e
) 
```




<hr>



### function eventheapdelete 

```C++
void eventheapdelete (
    struct event ** heap,
    int heapsize,
    int eventnum
) 
```




<hr>



### function eventheapify 

```C++
void eventheapify (
    struct event ** heap,
    int heapsize,
    int eventnum
) 
```




<hr>



### function eventheapinsert 

```C++
void eventheapinsert (
    struct event ** heap,
    int heapsize,
    struct event * newevent
) 
```




<hr>



### function exactinit 

```C++
void exactinit () 
```




<hr>



### function fast\_expansion\_sum\_zeroelim 

```C++
int fast_expansion_sum_zeroelim (
    int elen,
    REAL * e,
    int flen,
    REAL * f,
    REAL * h
) 
```




<hr>



### function findcircumcenter 

```C++
void findcircumcenter (
    struct mesh * m,
    struct behavior * b,
    vertex torg,
    vertex tdest,
    vertex tapex,
    vertex circumcenter,
    REAL * xi,
    REAL * eta,
    int offcenter
) 
```




<hr>



### function finddirection 

```C++
enum finddirectionresult finddirection (
    struct mesh * m,
    struct behavior * b,
    struct otri * searchtri,
    vertex searchpoint
) 
```




<hr>



### function flip 

```C++
void flip (
    struct mesh * m,
    struct behavior * b,
    struct otri * flipedge
) 
```




<hr>



### function formskeleton 

```C++
void formskeleton (
    struct mesh * m,
    struct behavior * b,
    int * segmentlist,
    int * segmentmarkerlist,
    int numberofsegments
) 
```




<hr>



### function frontlocate 

```C++
struct splaynode * frontlocate (
    struct mesh * m,
    struct splaynode * splayroot,
    struct otri * bottommost,
    vertex searchvertex,
    struct otri * searchtri,
    int * farright
) 
```




<hr>



### function getvertex 

```C++
vertex getvertex (
    struct mesh * m,
    struct behavior * b,
    int number
) 
```




<hr>



### function highorder 

```C++
void highorder (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function incircle 

```C++
REAL incircle (
    struct mesh * m,
    struct behavior * b,
    vertex pa,
    vertex pb,
    vertex pc,
    vertex pd
) 
```




<hr>



### function incircleadapt 

```C++
REAL incircleadapt (
    vertex pa,
    vertex pb,
    vertex pc,
    vertex pd,
    REAL permanent
) 
```




<hr>



### function incrementaldelaunay 

```C++
long incrementaldelaunay (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function infecthull 

```C++
void infecthull (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function initializetrisubpools 

```C++
void initializetrisubpools (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function initializevertexpool 

```C++
void initializevertexpool (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function insertsegment 

```C++
void insertsegment (
    struct mesh * m,
    struct behavior * b,
    vertex endpoint1,
    vertex endpoint2,
    int newmark
) 
```




<hr>



### function insertsubseg 

```C++
void insertsubseg (
    struct mesh * m,
    struct behavior * b,
    struct otri * tri,
    int subsegmark
) 
```




<hr>



### function insertvertex 

```C++
enum insertvertexresult insertvertex (
    struct mesh * m,
    struct behavior * b,
    vertex newvertex,
    struct otri * searchtri,
    struct osub * splitseg,
    int segmentflaws,
    int triflaws
) 
```




<hr>



### function internalerror 

```C++
void internalerror () 
```




<hr>



### function locate 

```C++
enum locateresult locate (
    struct mesh * m,
    struct behavior * b,
    vertex searchpoint,
    struct otri * searchtri
) 
```




<hr>



### function makesubseg 

```C++
void makesubseg (
    struct mesh * m,
    struct osub * newsubseg
) 
```




<hr>



### function maketriangle 

```C++
void maketriangle (
    struct mesh * m,
    struct behavior * b,
    struct otri * newotri
) 
```




<hr>



### function makevertexmap 

```C++
void makevertexmap (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function markhull 

```C++
void markhull (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function mergehulls 

```C++
void mergehulls (
    struct mesh * m,
    struct behavior * b,
    struct otri * farleft,
    struct otri * innerleft,
    struct otri * innerright,
    struct otri * farright,
    int axis
) 
```




<hr>



### function nonregular 

```C++
REAL nonregular (
    struct mesh * m,
    struct behavior * b,
    vertex pa,
    vertex pb,
    vertex pc,
    vertex pd
) 
```




<hr>



### function numbernodes 

```C++
void numbernodes (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function orient3d 

```C++
REAL orient3d (
    struct mesh * m,
    struct behavior * b,
    vertex pa,
    vertex pb,
    vertex pc,
    vertex pd,
    REAL aheight,
    REAL bheight,
    REAL cheight,
    REAL dheight
) 
```




<hr>



### function orient3dadapt 

```C++
REAL orient3dadapt (
    vertex pa,
    vertex pb,
    vertex pc,
    vertex pd,
    REAL aheight,
    REAL bheight,
    REAL cheight,
    REAL dheight,
    REAL permanent
) 
```




<hr>



### function parsecommandline 

```C++
void parsecommandline (
    int argc,
    char ** argv,
    struct behavior * b
) 
```




<hr>



### function plague 

```C++
void plague (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function poolalloc 

```C++
VOID * poolalloc (
    struct memorypool * pool
) 
```




<hr>



### function pooldealloc 

```C++
void pooldealloc (
    struct memorypool * pool,
    VOID * dyingitem
) 
```




<hr>



### function pooldeinit 

```C++
void pooldeinit (
    struct memorypool * pool
) 
```




<hr>



### function poolinit 

```C++
void poolinit (
    struct memorypool * pool,
    int bytecount,
    int itemcount,
    int firstitemcount,
    int alignment
) 
```




<hr>



### function poolrestart 

```C++
void poolrestart (
    struct memorypool * pool
) 
```




<hr>



### function poolzero 

```C++
void poolzero (
    struct memorypool * pool
) 
```




<hr>



### function preciselocate 

```C++
enum locateresult preciselocate (
    struct mesh * m,
    struct behavior * b,
    vertex searchpoint,
    struct otri * searchtri,
    int stopatsubsegment
) 
```




<hr>



### function precisionerror 

```C++
void precisionerror () 
```




<hr>



### function printsubseg 

```C++
void printsubseg (
    struct mesh * m,
    struct behavior * b,
    struct osub * s
) 
```




<hr>



### function printtriangle 

```C++
void printtriangle (
    struct mesh * m,
    struct behavior * b,
    struct otri * t
) 
```




<hr>



### function quality\_statistics 

```C++
void quality_statistics (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function randomnation 

```C++
unsigned long randomnation (
    unsigned int choices
) 
```




<hr>



### function reconstruct 

```C++
int reconstruct (
    struct mesh * m,
    struct behavior * b,
    int * trianglelist,
    REAL * triangleattriblist,
    REAL * trianglearealist,
    int elements,
    int corners,
    int attribs,
    int * segmentlist,
    int * segmentmarkerlist,
    int numberofsegments
) 
```




<hr>



### function regionplague 

```C++
void regionplague (
    struct mesh * m,
    struct behavior * b,
    REAL attribute,
    REAL area
) 
```




<hr>



### function removebox 

```C++
long removebox (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function removeghosts 

```C++
long removeghosts (
    struct mesh * m,
    struct behavior * b,
    struct otri * startghost
) 
```




<hr>



### function report 

```C++
void report (
    struct triangulateio * io,
    int markers,
    int reporttriangles,
    int reportneighbors,
    int reportsegments,
    int reportedges,
    int reportnorms
) 
```




<hr>



### function rightofhyperbola 

```C++
int rightofhyperbola (
    struct mesh * m,
    struct otri * fronttri,
    vertex newsite
) 
```




<hr>



### function scale\_expansion\_zeroelim 

```C++
int scale_expansion_zeroelim (
    int elen,
    REAL * e,
    REAL b,
    REAL * h
) 
```




<hr>



### function scoutsegment 

```C++
int scoutsegment (
    struct mesh * m,
    struct behavior * b,
    struct otri * searchtri,
    vertex endpoint2,
    int newmark
) 
```




<hr>



### function segmentintersection 

```C++
void segmentintersection (
    struct mesh * m,
    struct behavior * b,
    struct otri * splittri,
    struct osub * splitsubseg,
    vertex endpoint2
) 
```




<hr>



### function splay 

```C++
struct splaynode * splay (
    struct mesh * m,
    struct splaynode * splaytree,
    vertex searchpoint,
    struct otri * searchtri
) 
```




<hr>



### function splayinsert 

```C++
struct splaynode * splayinsert (
    struct mesh * m,
    struct splaynode * splayroot,
    struct otri * newkey,
    vertex searchpoint
) 
```




<hr>



### function splitencsegs 

```C++
void splitencsegs (
    struct mesh * m,
    struct behavior * b,
    int triflaws
) 
```




<hr>



### function splittriangle 

```C++
void splittriangle (
    struct mesh * m,
    struct behavior * b,
    struct badtriang * badtri
) 
```




<hr>



### function statistics 

```C++
void statistics (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function subsegdealloc 

```C++
void subsegdealloc (
    struct mesh * m,
    subseg * dyingsubseg
) 
```




<hr>



### function subsegtraverse 

```C++
subseg * subsegtraverse (
    struct mesh * m
) 
```




<hr>



### function sweeplinedelaunay 

```C++
long sweeplinedelaunay (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function tallyencs 

```C++
void tallyencs (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function tallyfaces 

```C++
void tallyfaces (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function testtriangle 

```C++
void testtriangle (
    struct mesh * m,
    struct behavior * b,
    struct otri * testtri
) 
```




<hr>



### function transfernodes 

```C++
void transfernodes (
    struct mesh * m,
    struct behavior * b,
    REAL * pointlist,
    REAL * pointattriblist,
    int * pointmarkerlist,
    int numberofpoints,
    int numberofpointattribs
) 
```




<hr>



### function traversalinit 

```C++
void traversalinit (
    struct memorypool * pool
) 
```




<hr>



### function traverse 

```C++
VOID * traverse (
    struct memorypool * pool
) 
```




<hr>



### function triangledealloc 

```C++
void triangledealloc (
    struct mesh * m,
    triangle * dyingtriangle
) 
```




<hr>



### function triangledeinit 

```C++
void triangledeinit (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function triangleinit 

```C++
void triangleinit (
    struct mesh * m
) 
```




<hr>



### function triangletraverse 

```C++
triangle * triangletraverse (
    struct mesh * m
) 
```




<hr>



### function triangulate 

```C++
void triangulate (
    char * triswitches,
    struct triangulateio * in,
    struct triangulateio * out,
    struct triangulateio * vorout
) 
```




<hr>



### function triangulatepolygon 

```C++
void triangulatepolygon (
    struct mesh * m,
    struct behavior * b,
    struct otri * firstedge,
    struct otri * lastedge,
    int edgecount,
    int doflip,
    int triflaws
) 
```




<hr>



### function triexit 

```C++
void triexit (
    int status
) 
```




<hr>



### function trifree 

```C++
void trifree (
    VOID * memptr
) 
```




<hr>



### function trimalloc 

```C++
VOID * trimalloc (
    int size
) 
```




<hr>



### function triunsuitable 

```C++
int triunsuitable (
    vertex triorg,
    vertex tridest,
    vertex triapex,
    REAL area
) 
```




<hr>



### function undovertex 

```C++
void undovertex (
    struct mesh * m,
    struct behavior * b
) 
```




<hr>



### function unflip 

```C++
void unflip (
    struct mesh * m,
    struct behavior * b,
    struct otri * flipedge
) 
```




<hr>



### function vertexdealloc 

```C++
void vertexdealloc (
    struct mesh * m,
    vertex dyingvertex
) 
```




<hr>



### function vertexmedian 

```C++
void vertexmedian (
    vertex * sortarray,
    int arraysize,
    int median,
    int axis
) 
```




<hr>



### function vertexsort 

```C++
void vertexsort (
    vertex * sortarray,
    int arraysize
) 
```




<hr>



### function vertextraverse 

```C++
vertex vertextraverse (
    struct mesh * m
) 
```




<hr>



### function writeedges 

```C++
void writeedges (
    struct mesh * m,
    struct behavior * b,
    int ** edgelist,
    int ** edgemarkerlist
) 
```




<hr>



### function writeelements 

```C++
void writeelements (
    struct mesh * m,
    struct behavior * b,
    int ** trianglelist,
    REAL ** triangleattriblist
) 
```




<hr>



### function writeneighbors 

```C++
void writeneighbors (
    struct mesh * m,
    struct behavior * b,
    int ** neighborlist
) 
```




<hr>



### function writenodes 

```C++
void writenodes (
    struct mesh * m,
    struct behavior * b,
    REAL ** pointlist,
    REAL ** pointattriblist,
    int ** pointmarkerlist
) 
```




<hr>



### function writepoly 

```C++
void writepoly (
    struct mesh * m,
    struct behavior * b,
    int ** segmentlist,
    int ** segmentmarkerlist
) 
```




<hr>



### function writevoronoi 

```C++
void writevoronoi (
    struct mesh * m,
    struct behavior * b,
    REAL ** vpointlist,
    REAL ** vpointattriblist,
    int ** vpointmarkerlist,
    int ** vedgelist,
    int ** vedgemarkerlist,
    REAL ** vnormlist
) 
```




<hr>
## Macro Definition Documentation





### define Absolute 

```C++
#define Absolute (
    a
) `((a) >= 0.0 ? (a) : -(a))`
```




<hr>



### define BADSUBSEGPERBLOCK 

```C++
#define BADSUBSEGPERBLOCK `252`
```




<hr>



### define BADTRIPERBLOCK 

```C++
#define BADTRIPERBLOCK `4092`
```




<hr>



### define DEADVERTEX 

```C++
#define DEADVERTEX `-32768`
```




<hr>



### define FILENAMESIZE 

```C++
#define FILENAMESIZE `2048`
```




<hr>



### define FLIPSTACKERPERBLOCK 

```C++
#define FLIPSTACKERPERBLOCK `252`
```




<hr>



### define FREEVERTEX 

```C++
#define FREEVERTEX `2`
```




<hr>



### define Fast\_Two\_Sum 

```C++
#define Fast_Two_Sum (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Fast\_Two\_Sum\_Tail 

```C++
#define Fast_Two_Sum_Tail (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define INEXACT 

```C++
#define INEXACT `/* Nothing */`
```




<hr>



### define INPUTLINESIZE 

```C++
#define INPUTLINESIZE `1024`
```




<hr>



### define INPUTVERTEX 

```C++
#define INPUTVERTEX `0`
```




<hr>



### define ONETHIRD 

```C++
#define ONETHIRD `0.333333333333333333333333333333333333333333333333333333333333`
```




<hr>



### define PI 

```C++
#define PI `3.141592653589793238462643383279502884197169399375105820974944592308`
```




<hr>



### define REAL 

```C++
#define REAL `double`
```




<hr>



### define SAMPLEFACTOR 

```C++
#define SAMPLEFACTOR `11`
```




<hr>



### define SAMPLERATE 

```C++
#define SAMPLERATE `10`
```




<hr>



### define SEGMENTVERTEX 

```C++
#define SEGMENTVERTEX `1`
```




<hr>



### define SPLAYNODEPERBLOCK 

```C++
#define SPLAYNODEPERBLOCK `508`
```




<hr>



### define SQUAREROOTTWO 

```C++
#define SQUAREROOTTWO `1.4142135623730950488016887242096980785696718753769480732`
```




<hr>



### define STARTINDEX 

```C++
#define STARTINDEX `0`
```




<hr>



### define SUBSEGPERBLOCK 

```C++
#define SUBSEGPERBLOCK `508       /* Number of subsegments allocated at once. */`
```




<hr>



### define Split 

```C++
#define Split (
    a,
    ahi,
    alo
) `/* multi line expression */`
```




<hr>



### define Square 

```C++
#define Square (
    a,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Square\_Tail 

```C++
#define Square_Tail (
    a,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define TRILIBRARY 

```C++
#define TRILIBRARY `1`
```




<hr>



### define TRIPERBLOCK 

```C++
#define TRIPERBLOCK `4092           /* Number of triangles allocated at once. */`
```




<hr>



### define Two\_Diff 

```C++
#define Two_Diff (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Diff\_Tail 

```C++
#define Two_Diff_Tail (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_One\_Diff 

```C++
#define Two_One_Diff (
    a1,
    a0,
    b,
    x2,
    x1,
    x0
) `/* multi line expression */`
```




<hr>



### define Two\_One\_Product 

```C++
#define Two_One_Product (
    a1,
    a0,
    b,
    x3,
    x2,
    x1,
    x0
) `/* multi line expression */`
```




<hr>



### define Two\_One\_Sum 

```C++
#define Two_One_Sum (
    a1,
    a0,
    b,
    x2,
    x1,
    x0
) `/* multi line expression */`
```




<hr>



### define Two\_Product 

```C++
#define Two_Product (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Product\_Presplit 

```C++
#define Two_Product_Presplit (
    a,
    b,
    bhi,
    blo,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Product\_Tail 

```C++
#define Two_Product_Tail (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Sum 

```C++
#define Two_Sum (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Sum\_Tail 

```C++
#define Two_Sum_Tail (
    a,
    b,
    x,
    y
) `/* multi line expression */`
```




<hr>



### define Two\_Two\_Diff 

```C++
#define Two_Two_Diff (
    a1,
    a0,
    b1,
    b0,
    x3,
    x2,
    x1,
    x0
) `/* multi line expression */`
```




<hr>



### define Two\_Two\_Sum 

```C++
#define Two_Two_Sum (
    a1,
    a0,
    b1,
    b0,
    x3,
    x2,
    x1,
    x0
) `/* multi line expression */`
```




<hr>



### define UNDEADVERTEX 

```C++
#define UNDEADVERTEX `-32767`
```




<hr>



### define VERTEXPERBLOCK 

```C++
#define VERTEXPERBLOCK `4092         /* Number of vertices allocated at once. */`
```




<hr>



### define VIRUSPERBLOCK 

```C++
#define VIRUSPERBLOCK `1020   /* Number of virus triangles allocated at once. */`
```




<hr>



### define VOID 

```C++
#define VOID `int`
```




<hr>



### define apex 

```C++
#define apex (
    otri,
    vertexptr
) `vertexptr = (vertex) ( otri ).tri[( otri ).orient + 3]`
```




<hr>



### define areabound 

```C++
#define areabound (
    otri
) `((REAL *) ( otri ).tri)[m->areaboundindex]`
```




<hr>



### define bond 

```C++
#define bond (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define deadsubseg 

```C++
#define deadsubseg (
    sub
) `((sub)[1] == (subseg) NULL)`
```




<hr>



### define deadtri 

```C++
#define deadtri (
    tria
) `((tria)[1] == (triangle) NULL)`
```




<hr>



### define decode 

```C++
#define decode (
    ptr,
    otri
) `/* multi line expression */`
```




<hr>



### define dest 

```C++
#define dest (
    otri,
    vertexptr
) `vertexptr = (vertex) ( otri ).tri[minus1mod3[( otri ).orient] + 3]`
```




<hr>



### define dissolve 

```C++
#define dissolve (
    otri
) `( otri ).tri[( otri ).orient] = (triangle) m->dummytri`
```




<hr>



### define dnext 

```C++
#define dnext (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define dnextself 

```C++
#define dnextself (
    otri
) `/* multi line expression */`
```




<hr>



### define dprev 

```C++
#define dprev (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define dprevself 

```C++
#define dprevself (
    otri
) `/* multi line expression */`
```




<hr>



### define elemattribute 

```C++
#define elemattribute (
    otri,
    attnum
) `((REAL *) ( otri ).tri)[m->elemattribindex + (attnum)]`
```




<hr>



### define encode 

```C++
#define encode (
    otri
) `(triangle) ((unsigned long) ( otri ).tri | (unsigned long) ( otri ).orient)`
```




<hr>



### define infect 

```C++
#define infect (
    otri
) `/* multi line expression */`
```




<hr>



### define infected 

```C++
#define infected (
    otri
) `(((unsigned long) ( otri ).tri[6] & (unsigned long) 2l) != 0l)`
```




<hr>



### define killsubseg 

```C++
#define killsubseg (
    sub
) `/* multi line expression */`
```




<hr>



### define killtri 

```C++
#define killtri (
    tria
) `/* multi line expression */`
```




<hr>



### define lnext 

```C++
#define lnext (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define lnextself 

```C++
#define lnextself (
    otri
) `( otri ).orient = plus1mod3[( otri ).orient]`
```




<hr>



### define lprev 

```C++
#define lprev (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define lprevself 

```C++
#define lprevself (
    otri
) `( otri ).orient = minus1mod3[( otri ).orient]`
```




<hr>



### define mark 

```C++
#define mark (
    osub
) `(* (int *) (( osub ).ss + 8))`
```




<hr>



### define onext 

```C++
#define onext (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define onextself 

```C++
#define onextself (
    otri
) `/* multi line expression */`
```




<hr>



### define oprev 

```C++
#define oprev (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define oprevself 

```C++
#define oprevself (
    otri
) `/* multi line expression */`
```




<hr>



### define org 

```C++
#define org (
    otri,
    vertexptr
) `vertexptr = (vertex) ( otri ).tri[plus1mod3[( otri ).orient] + 3]`
```




<hr>



### define otricopy 

```C++
#define otricopy (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define otriequal 

```C++
#define otriequal (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define rnext 

```C++
#define rnext (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define rnextself 

```C++
#define rnextself (
    otri
) `/* multi line expression */`
```




<hr>



### define rprev 

```C++
#define rprev (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define rprevself 

```C++
#define rprevself (
    otri
) `/* multi line expression */`
```




<hr>



### define sbond 

```C++
#define sbond (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define sdecode 

```C++
#define sdecode (
    sptr,
    osub
) `/* multi line expression */`
```




<hr>



### define sdest 

```C++
#define sdest (
    osub,
    vertexptr
) `vertexptr = (vertex) ( osub ).ss[3 - ( osub ).ssorient]`
```




<hr>



### define sdissolve 

```C++
#define sdissolve (
    osub
) `( osub ).ss[( osub ).ssorient] = (subseg) m->dummysub`
```




<hr>



### define segdest 

```C++
#define segdest (
    osub,
    vertexptr
) `vertexptr = (vertex) ( osub ).ss[5 - ( osub ).ssorient]`
```




<hr>



### define segorg 

```C++
#define segorg (
    osub,
    vertexptr
) `vertexptr = (vertex) ( osub ).ss[4 + ( osub ).ssorient]`
```




<hr>



### define sencode 

```C++
#define sencode (
    osub
) `(subseg) ((unsigned long) ( osub ).ss | (unsigned long) ( osub ).ssorient)`
```




<hr>



### define setapex 

```C++
#define setapex (
    otri,
    vertexptr
) `( otri ).tri[( otri ).orient + 3] = (triangle) vertexptr`
```




<hr>



### define setareabound 

```C++
#define setareabound (
    otri,
    value
) `((REAL *) ( otri ).tri)[m->areaboundindex] = value`
```




<hr>



### define setdest 

```C++
#define setdest (
    otri,
    vertexptr
) `( otri ).tri[minus1mod3[( otri ).orient] + 3] = (triangle) vertexptr`
```




<hr>



### define setelemattribute 

```C++
#define setelemattribute (
    otri,
    attnum,
    value
) `((REAL *) ( otri ).tri)[m->elemattribindex + (attnum)] = value`
```




<hr>



### define setmark 

```C++
#define setmark (
    osub,
    value
) `* (int *) (( osub ).ss + 8) = value`
```




<hr>



### define setorg 

```C++
#define setorg (
    otri,
    vertexptr
) `( otri ).tri[plus1mod3[( otri ).orient] + 3] = (triangle) vertexptr`
```




<hr>



### define setsdest 

```C++
#define setsdest (
    osub,
    vertexptr
) `( osub ).ss[3 - ( osub ).ssorient] = (subseg) vertexptr`
```




<hr>



### define setsegdest 

```C++
#define setsegdest (
    osub,
    vertexptr
) `( osub ).ss[5 - ( osub ).ssorient] = (subseg) vertexptr`
```




<hr>



### define setsegorg 

```C++
#define setsegorg (
    osub,
    vertexptr
) `( osub ).ss[4 + ( osub ).ssorient] = (subseg) vertexptr`
```




<hr>



### define setsorg 

```C++
#define setsorg (
    osub,
    vertexptr
) `( osub ).ss[2 + ( osub ).ssorient] = (subseg) vertexptr`
```




<hr>



### define setvertex2tri 

```C++
#define setvertex2tri (
    vx,
    value
) `((triangle *) (vx))[m->vertex2triindex] = value`
```




<hr>



### define setvertexmark 

```C++
#define setvertexmark (
    vx,
    value
) `((int *) (vx))[m->vertexmarkindex] = value`
```




<hr>



### define setvertextype 

```C++
#define setvertextype (
    vx,
    value
) `((int *) (vx))[m->vertexmarkindex + 1] = value`
```




<hr>



### define snext 

```C++
#define snext (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define snextself 

```C++
#define snextself (
    osub
) `/* multi line expression */`
```




<hr>



### define sorg 

```C++
#define sorg (
    osub,
    vertexptr
) `vertexptr = (vertex) ( osub ).ss[2 + ( osub ).ssorient]`
```




<hr>



### define spivot 

```C++
#define spivot (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define spivotself 

```C++
#define spivotself (
    osub
) `/* multi line expression */`
```




<hr>



### define ssym 

```C++
#define ssym (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define ssymself 

```C++
#define ssymself (
    osub
) `( osub ).ssorient = 1 - ( osub ).ssorient`
```




<hr>



### define stdissolve 

```C++
#define stdissolve (
    osub
) `( osub ).ss[6 + ( osub ).ssorient] = (subseg) m->dummytri`
```




<hr>



### define stpivot 

```C++
#define stpivot (
    osub,
    otri
) `/* multi line expression */`
```




<hr>



### define subsegcopy 

```C++
#define subsegcopy (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define subsegequal 

```C++
#define subsegequal (
    osub1,
    osub2
) `/* multi line expression */`
```




<hr>



### define sym 

```C++
#define sym (
    otri1,
    otri2
) `/* multi line expression */`
```




<hr>



### define symself 

```C++
#define symself (
    otri
) `/* multi line expression */`
```




<hr>



### define tsbond 

```C++
#define tsbond (
    otri,
    osub
) `/* multi line expression */`
```




<hr>



### define tsdissolve 

```C++
#define tsdissolve (
    otri
) `( otri ).tri[6 + ( otri ).orient] = (triangle) m->dummysub`
```




<hr>



### define tspivot 

```C++
#define tspivot (
    otri,
    osub
) `/* multi line expression */`
```




<hr>



### define uninfect 

```C++
#define uninfect (
    otri
) `/* multi line expression */`
```




<hr>



### define vertex2tri 

```C++
#define vertex2tri (
    vx
) `((triangle *) (vx))[m->vertex2triindex]`
```




<hr>



### define vertexmark 

```C++
#define vertexmark (
    vx
) `((int *) (vx))[m->vertexmarkindex]`
```




<hr>



### define vertextype 

```C++
#define vertextype (
    vx
) `((int *) (vx))[m->vertexmarkindex + 1]`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/triangle.cpp`

