

# File MultiProgressBar.h



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**MultiProgressBar.h**](MultiProgressBar_8h.md)

[Go to the source code of this file](MultiProgressBar_8h_source.md)

_Definition of progress bar._ [More...](#detailed-description)

* `#include <iostream>`
* `#include <vector>`
* `#include <cmath>`
* `#include <string>`
* `#include <sys/ioctl.h>`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**MultiProgressBar**](classMultiProgressBar.md) <br> |

















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**CLEAR**](MultiProgressBar_8h.md#define-clear) () `printf("\033[2J")`<br> |
| define  | [**COLOR\_BAR\_BLUE**](MultiProgressBar_8h.md#define-color_bar_blue)  `1`<br> |
| define  | [**COLOR\_BAR\_GREEN**](MultiProgressBar_8h.md#define-color_bar_green)  `2`<br> |
| define  | [**COLOR\_BAR\_PURPLE**](MultiProgressBar_8h.md#define-color_bar_purple)  `0`<br> |
| define  | [**COLOR\_BAR\_RED**](MultiProgressBar_8h.md#define-color_bar_red)  `4`<br> |
| define  | [**COLOR\_BAR\_YELLOW**](MultiProgressBar_8h.md#define-color_bar_yellow)  `3`<br> |
| define  | [**HIDE\_CURSOR**](MultiProgressBar_8h.md#define-hide_cursor) () `printf("\033[?25l")`<br> |
| define  | [**HIGHT\_LIGHT**](MultiProgressBar_8h.md#define-hight_light) () `printf("\033[7m")`<br> |
| define  | [**MOVEDOWN**](MultiProgressBar_8h.md#define-movedown) (x) `printf("\033[%dB", (x))`<br> |
| define  | [**MOVELEFT**](MultiProgressBar_8h.md#define-moveleft) (y) `printf("\033[%dD", (y))`<br> |
| define  | [**MOVERIGHT**](MultiProgressBar_8h.md#define-moveright) (y) `printf("\033[%dC",(y))`<br> |
| define  | [**MOVETO**](MultiProgressBar_8h.md#define-moveto) (x, y) `printf("\033[%d;%dH", (x), (y))`<br> |
| define  | [**MOVEUP**](MultiProgressBar_8h.md#define-moveup) (x) `printf("\033[%dA", (x))`<br> |
| define  | [**MOVEUP**](MultiProgressBar_8h.md#define-moveup) (x) `printf("\033[%dA", (x))`<br> |
| define  | [**RESET\_CURSOR**](MultiProgressBar_8h.md#define-reset_cursor) () `printf("\033[H")`<br> |
| define  | [**SHOW\_CURSOR**](MultiProgressBar_8h.md#define-show_cursor) () `printf("\033[?25h")`<br> |
| define  | [**UN\_HIGHT\_LIGHT**](MultiProgressBar_8h.md#define-un_hight_light) () `printf("\033[27m")`<br> |

## Detailed Description




**Author:**

Zhikui Guo ([zhikuiguo@live.cn](mailto:zhikuiguo@live.cn)) 




**Version:**

1.0 




**Date:**

2019-09-03




**Copyright:**

Copyright (c) 2019 





    
## Macro Definition Documentation





### define CLEAR 

```C++
#define CLEAR (
    
) `printf("\033[2J")`
```




<hr>



### define COLOR\_BAR\_BLUE 

```C++
#define COLOR_BAR_BLUE `1`
```




<hr>



### define COLOR\_BAR\_GREEN 

```C++
#define COLOR_BAR_GREEN `2`
```




<hr>



### define COLOR\_BAR\_PURPLE 

```C++
#define COLOR_BAR_PURPLE `0`
```




<hr>



### define COLOR\_BAR\_RED 

```C++
#define COLOR_BAR_RED `4`
```




<hr>



### define COLOR\_BAR\_YELLOW 

```C++
#define COLOR_BAR_YELLOW `3`
```




<hr>



### define HIDE\_CURSOR 

```C++
#define HIDE_CURSOR (
    
) `printf("\033[?25l")`
```




<hr>



### define HIGHT\_LIGHT 

```C++
#define HIGHT_LIGHT (
    
) `printf("\033[7m")`
```




<hr>



### define MOVEDOWN 

```C++
#define MOVEDOWN (
    x
) `printf("\033[%dB", (x))`
```




<hr>



### define MOVELEFT 

```C++
#define MOVELEFT (
    y
) `printf("\033[%dD", (y))`
```




<hr>



### define MOVERIGHT 

```C++
#define MOVERIGHT (
    y
) `printf("\033[%dC",(y))`
```




<hr>



### define MOVETO 

```C++
#define MOVETO (
    x,
    y
) `printf("\033[%d;%dH", (x), (y))`
```




<hr>



### define MOVEUP 

```C++
#define MOVEUP (
    x
) `printf("\033[%dA", (x))`
```




<hr>



### define MOVEUP 

```C++
#define MOVEUP (
    x
) `printf("\033[%dA", (x))`
```




<hr>



### define RESET\_CURSOR 

```C++
#define RESET_CURSOR (
    
) `printf("\033[H")`
```




<hr>



### define SHOW\_CURSOR 

```C++
#define SHOW_CURSOR (
    
) `printf("\033[?25h")`
```




<hr>



### define UN\_HIGHT\_LIGHT 

```C++
#define UN_HIGHT_LIGHT (
    
) `printf("\033[27m")`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/MultiProgressBar.h`

