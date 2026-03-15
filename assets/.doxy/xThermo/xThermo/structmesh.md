

# Struct mesh



[**ClassList**](annotated.md) **>** [**mesh**](structmesh.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  int | [**areaboundindex**](#variable-areaboundindex)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**badsubsegs**](#variable-badsubsegs)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**badtriangles**](#variable-badtriangles)  <br> |
|  int | [**checkquality**](#variable-checkquality)  <br> |
|  int | [**checksegments**](#variable-checksegments)  <br> |
|  long | [**circletopcount**](#variable-circletopcount)  <br> |
|  long | [**circumcentercount**](#variable-circumcentercount)  <br> |
|  long | [**counterclockcount**](#variable-counterclockcount)  <br> |
|  subseg \* | [**dummysub**](#variable-dummysub)  <br> |
|  subseg \* | [**dummysubbase**](#variable-dummysubbase)  <br> |
|  triangle \* | [**dummytri**](#variable-dummytri)  <br> |
|  triangle \* | [**dummytribase**](#variable-dummytribase)  <br> |
|  long | [**edges**](#variable-edges)  <br> |
|  int | [**eextras**](#variable-eextras)  <br> |
|  int | [**elemattribindex**](#variable-elemattribindex)  <br> |
|  int | [**firstnonemptyq**](#variable-firstnonemptyq)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**flipstackers**](#variable-flipstackers)  <br> |
|  int | [**highorderindex**](#variable-highorderindex)  <br> |
|  int | [**holes**](#variable-holes)  <br> |
|  long | [**hullsize**](#variable-hullsize)  <br> |
|  long | [**hyperbolacount**](#variable-hyperbolacount)  <br> |
|  long | [**incirclecount**](#variable-incirclecount)  <br> |
|  int | [**inelements**](#variable-inelements)  <br> |
|  vertex | [**infvertex1**](#variable-infvertex1)  <br> |
|  vertex | [**infvertex2**](#variable-infvertex2)  <br> |
|  vertex | [**infvertex3**](#variable-infvertex3)  <br> |
|  int | [**insegments**](#variable-insegments)  <br> |
|  int | [**invertices**](#variable-invertices)  <br> |
|  struct [**flipstacker**](structflipstacker.md) \* | [**lastflip**](#variable-lastflip)  <br> |
|  int | [**mesh\_dim**](#variable-mesh_dim)  <br> |
|  int | [**nextnonemptyq**](#variable-nextnonemptyq)  <br> |
|  int | [**nextras**](#variable-nextras)  <br> |
|  long | [**orient3dcount**](#variable-orient3dcount)  <br> |
|  struct [**badtriang**](structbadtriang.md) \* | [**queuefront**](#variable-queuefront)  <br> |
|  struct [**badtriang**](structbadtriang.md) \* | [**queuetail**](#variable-queuetail)  <br> |
|  int | [**readnodefile**](#variable-readnodefile)  <br> |
|  struct [**otri**](structotri.md) | [**recenttri**](#variable-recenttri)  <br> |
|  int | [**regions**](#variable-regions)  <br> |
|  long | [**samples**](#variable-samples)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**splaynodes**](#variable-splaynodes)  <br> |
|  int | [**steinerleft**](#variable-steinerleft)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**subsegs**](#variable-subsegs)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**triangles**](#variable-triangles)  <br> |
|  int | [**undeads**](#variable-undeads)  <br> |
|  int | [**vertex2triindex**](#variable-vertex2triindex)  <br> |
|  int | [**vertexmarkindex**](#variable-vertexmarkindex)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**vertices**](#variable-vertices)  <br> |
|  struct [**memorypool**](structmemorypool.md) | [**viri**](#variable-viri)  <br> |
|  REAL | [**xmax**](#variable-xmax)  <br> |
|  REAL | [**xmin**](#variable-xmin)  <br> |
|  REAL | [**xminextreme**](#variable-xminextreme)  <br> |
|  REAL | [**ymax**](#variable-ymax)  <br> |
|  REAL | [**ymin**](#variable-ymin)  <br> |












































## Public Attributes Documentation




### variable areaboundindex 

```C++
int mesh::areaboundindex;
```




<hr>



### variable badsubsegs 

```C++
struct memorypool mesh::badsubsegs;
```




<hr>



### variable badtriangles 

```C++
struct memorypool mesh::badtriangles;
```




<hr>



### variable checkquality 

```C++
int mesh::checkquality;
```




<hr>



### variable checksegments 

```C++
int mesh::checksegments;
```




<hr>



### variable circletopcount 

```C++
long mesh::circletopcount;
```




<hr>



### variable circumcentercount 

```C++
long mesh::circumcentercount;
```




<hr>



### variable counterclockcount 

```C++
long mesh::counterclockcount;
```




<hr>



### variable dummysub 

```C++
subseg* mesh::dummysub;
```




<hr>



### variable dummysubbase 

```C++
subseg* mesh::dummysubbase;
```




<hr>



### variable dummytri 

```C++
triangle* mesh::dummytri;
```




<hr>



### variable dummytribase 

```C++
triangle* mesh::dummytribase;
```




<hr>



### variable edges 

```C++
long mesh::edges;
```




<hr>



### variable eextras 

```C++
int mesh::eextras;
```




<hr>



### variable elemattribindex 

```C++
int mesh::elemattribindex;
```




<hr>



### variable firstnonemptyq 

```C++
int mesh::firstnonemptyq;
```




<hr>



### variable flipstackers 

```C++
struct memorypool mesh::flipstackers;
```




<hr>



### variable highorderindex 

```C++
int mesh::highorderindex;
```




<hr>



### variable holes 

```C++
int mesh::holes;
```




<hr>



### variable hullsize 

```C++
long mesh::hullsize;
```




<hr>



### variable hyperbolacount 

```C++
long mesh::hyperbolacount;
```




<hr>



### variable incirclecount 

```C++
long mesh::incirclecount;
```




<hr>



### variable inelements 

```C++
int mesh::inelements;
```




<hr>



### variable infvertex1 

```C++
vertex mesh::infvertex1;
```




<hr>



### variable infvertex2 

```C++
vertex mesh::infvertex2;
```




<hr>



### variable infvertex3 

```C++
vertex mesh::infvertex3;
```




<hr>



### variable insegments 

```C++
int mesh::insegments;
```




<hr>



### variable invertices 

```C++
int mesh::invertices;
```




<hr>



### variable lastflip 

```C++
struct flipstacker* mesh::lastflip;
```




<hr>



### variable mesh\_dim 

```C++
int mesh::mesh_dim;
```




<hr>



### variable nextnonemptyq 

```C++
int mesh::nextnonemptyq[4096];
```




<hr>



### variable nextras 

```C++
int mesh::nextras;
```




<hr>



### variable orient3dcount 

```C++
long mesh::orient3dcount;
```




<hr>



### variable queuefront 

```C++
struct badtriang* mesh::queuefront[4096];
```




<hr>



### variable queuetail 

```C++
struct badtriang* mesh::queuetail[4096];
```




<hr>



### variable readnodefile 

```C++
int mesh::readnodefile;
```




<hr>



### variable recenttri 

```C++
struct otri mesh::recenttri;
```




<hr>



### variable regions 

```C++
int mesh::regions;
```




<hr>



### variable samples 

```C++
long mesh::samples;
```




<hr>



### variable splaynodes 

```C++
struct memorypool mesh::splaynodes;
```




<hr>



### variable steinerleft 

```C++
int mesh::steinerleft;
```




<hr>



### variable subsegs 

```C++
struct memorypool mesh::subsegs;
```




<hr>



### variable triangles 

```C++
struct memorypool mesh::triangles;
```




<hr>



### variable undeads 

```C++
int mesh::undeads;
```




<hr>



### variable vertex2triindex 

```C++
int mesh::vertex2triindex;
```




<hr>



### variable vertexmarkindex 

```C++
int mesh::vertexmarkindex;
```




<hr>



### variable vertices 

```C++
struct memorypool mesh::vertices;
```




<hr>



### variable viri 

```C++
struct memorypool mesh::viri;
```




<hr>



### variable xmax 

```C++
REAL mesh::xmax;
```




<hr>



### variable xmin 

```C++
REAL mesh::xmin;
```




<hr>



### variable xminextreme 

```C++
REAL mesh::xminextreme;
```




<hr>



### variable ymax 

```C++
REAL mesh::ymax;
```




<hr>



### variable ymin 

```C++
REAL mesh::ymin;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/triangle.cpp`

