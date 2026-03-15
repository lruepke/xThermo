

# File stdfunc.h



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**stdfunc.h**](test_2test__templateClass_2stdfunc_8h.md)

[Go to the source code of this file](test_2test__templateClass_2stdfunc_8h_source.md)



* `#include "stdio.h"`
* `#include <iostream>`
* `#include <map>`
* `#include <sys/stat.h>`
* `#include <cstring>`
* `#include <sstream>`
* `#include <iomanip>`
* `#include <cassert>`
































































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**ASSERT**](test_2test__templateClass_2stdfunc_8h.md#define-assert) (expression, info) `{if(!(expression))std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_RED&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl; assert(expression);}`<br> |
| define  | [**COLOR\_BLUE**](test_2test__templateClass_2stdfunc_8h.md#define-color_blue)  `"\033[34m"`<br> |
| define  | [**COLOR\_DEFAULT**](test_2test__templateClass_2stdfunc_8h.md#define-color_default)  `"\033[0m"`<br> |
| define  | [**COLOR\_GREEN**](test_2test__templateClass_2stdfunc_8h.md#define-color_green)  `"\033[32m"`<br> |
| define  | [**COLOR\_PURPLE**](test_2test__templateClass_2stdfunc_8h.md#define-color_purple)  `"\033[35m"`<br> |
| define  | [**COLOR\_RED**](test_2test__templateClass_2stdfunc_8h.md#define-color_red)  `"\033[31m"`<br> |
| define  | [**COLOR\_YELLOW**](test_2test__templateClass_2stdfunc_8h.md#define-color_yellow)  `"\033[33m"`<br> |
| define  | [**ERROR**](test_2test__templateClass_2stdfunc_8h.md#define-error) (info) `{std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_RED&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl; exit(0);}`<br> |
| define  | [**ERROR\_COUT**](test_2test__templateClass_2stdfunc_8h.md#define-error_cout)  `"["&lt;&lt;"\033[31mError: "&lt;&lt;"\033[0m] "`<br> |
| define  | [**STATUS**](test_2test__templateClass_2stdfunc_8h.md#define-status) (info) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_color**](test_2test__templateClass_2stdfunc_8h.md#define-status_color) (info, color) `std::cout&lt;&lt;"--  "&lt;&lt;color&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_system\_time**](test_2test__templateClass_2stdfunc_8h.md#define-status_system_time) (info, duration) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;", time: "&lt;&lt; double(duration) \* std::chrono::microseconds::period::num / std::chrono::microseconds::period::den&lt;&lt;" s"&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_time**](test_2test__templateClass_2stdfunc_8h.md#define-status_time) (info, time\_taken) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;", time: "&lt;&lt;(double)(time\_taken)/CLOCKS\_PER\_SEC&lt;&lt;" s"&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**WAIT**](test_2test__templateClass_2stdfunc_8h.md#define-wait) (where) `{std::cout&lt;&lt;"Waiting "&lt;&lt;where&lt;&lt; ". Enter to continue..." &lt;&lt; std::endl; std::string dummy; std::getline(std::cin, dummy);}`<br> |
| define  | [**WARNING**](test_2test__templateClass_2stdfunc_8h.md#define-warning) (info) `{std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_YELLOW&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;}`<br> |
| define  | [**WARN\_COUT**](test_2test__templateClass_2stdfunc_8h.md#define-warn_cout)  `"["&lt;&lt;"\033[33mWarning: "&lt;&lt;"\033[0m] "`<br> |

## Macro Definition Documentation





### define ASSERT 

```C++
#define ASSERT (
    expression,
    info
) `{if(!(expression))std::cout<<"--  "<<COLOR_RED<<info<<COLOR_DEFAULT<<std::endl; assert(expression);}`
```




<hr>



### define COLOR\_BLUE 

```C++
#define COLOR_BLUE `"\033[34m"`
```




<hr>



### define COLOR\_DEFAULT 

```C++
#define COLOR_DEFAULT `"\033[0m"`
```




<hr>



### define COLOR\_GREEN 

```C++
#define COLOR_GREEN `"\033[32m"`
```




<hr>



### define COLOR\_PURPLE 

```C++
#define COLOR_PURPLE `"\033[35m"`
```




<hr>



### define COLOR\_RED 

```C++
#define COLOR_RED `"\033[31m"`
```




<hr>



### define COLOR\_YELLOW 

```C++
#define COLOR_YELLOW `"\033[33m"`
```




<hr>



### define ERROR 

```C++
#define ERROR (
    info
) `{std::cout<<"--  "<<COLOR_RED<<info<<COLOR_DEFAULT<<std::endl; exit(0);}`
```




<hr>



### define ERROR\_COUT 

```C++
#define ERROR_COUT `"["<<"\033[31mError: "<<"\033[0m] "`
```




<hr>



### define STATUS 

```C++
#define STATUS (
    info
) `std::cout<<"--  "<<COLOR_GREEN<<info<<COLOR_DEFAULT<<std::endl;`
```




<hr>



### define STATUS\_color 

```C++
#define STATUS_color (
    info,
    color
) `std::cout<<"--  "<<color<<info<<COLOR_DEFAULT<<std::endl;`
```




<hr>



### define STATUS\_system\_time 

```C++
#define STATUS_system_time (
    info,
    duration
) `std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<< double(duration) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den<<" s"<<COLOR_DEFAULT<<std::endl;`
```




<hr>



### define STATUS\_time 

```C++
#define STATUS_time (
    info,
    time_taken
) `std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<<(double)(time_taken)/CLOCKS_PER_SEC<<" s"<<COLOR_DEFAULT<<std::endl;`
```




<hr>



### define WAIT 

```C++
#define WAIT (
    where
) `{std::cout<<"Waiting "<<where<< ". Enter to continue..." << std::endl; std::string dummy; std::getline(std::cin, dummy);}`
```




<hr>



### define WARNING 

```C++
#define WARNING (
    info
) `{std::cout<<"--  "<<COLOR_YELLOW<<info<<COLOR_DEFAULT<<std::endl;}`
```




<hr>



### define WARN\_COUT 

```C++
#define WARN_COUT `"["<<"\033[33mWarning: "<<"\033[0m] "`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/test/test_templateClass/stdfunc.h`

