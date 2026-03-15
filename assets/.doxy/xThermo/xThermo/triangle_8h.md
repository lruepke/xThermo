

# File triangle.h



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**triangle.h**](triangle_8h.md)

[Go to the source code of this file](triangle_8h_source.md)


















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**triangulateio**](structtriangulateio.md) <br> |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**report**](#function-report) (struct [**triangulateio**](structtriangulateio.md) \* io, int markers, int reporttriangles, int reportneighbors, int reportsegments, int reportedges, int reportnorms) <br> |
|  void | [**triangulate**](#function-triangulate) (char \* triswitches, struct [**triangulateio**](structtriangulateio.md) \* in, struct [**triangulateio**](structtriangulateio.md) \* out, struct [**triangulateio**](structtriangulateio.md) \* vorout) <br> |
|  void | [**trifree**](#function-trifree) (VOID \* memptr) <br> |



























## Macros

| Type | Name |
| ---: | :--- |
| define  | [**ANSI\_DECLARATORS**](triangle_8h.md#define-ansi_declarators)  <br> |
| define  | [**REAL**](triangle_8h.md#define-real)  `double`<br> |
| define  | [**VOID**](triangle_8h.md#define-void)  `int`<br> |

## Public Functions Documentation




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



### function trifree 

```C++
void trifree (
    VOID * memptr
) 
```




<hr>
## Macro Definition Documentation





### define ANSI\_DECLARATORS 

```C++
#define ANSI_DECLARATORS 
```




<hr>



### define REAL 

```C++
#define REAL `double`
```




<hr>



### define VOID 

```C++
#define VOID `int`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/triangle.h`

