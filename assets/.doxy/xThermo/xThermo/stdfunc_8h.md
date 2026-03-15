

# File stdfunc.h



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**stdfunc.h**](stdfunc_8h.md)

[Go to the source code of this file](stdfunc_8h_source.md)



* `#include <cstdio>`
* `#include <iostream>`
* `#include <map>`
* `#include <vector>`
* `#include <sys/stat.h>`
* `#include <cstring>`
* `#include <sstream>`
* `#include <iomanip>`
* `#include <cassert>`
* `#include <unistd.h>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |



















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**ASSERT**](stdfunc_8h.md#define-assert) (expression, info) `{if(!(expression))std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_RED&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl; assert(expression);}`<br> |
| define  | [**COLOR\_BLUE**](stdfunc_8h.md#define-color_blue)  `(isatty(STDOUT\_FILENO) == true ? "\033[34m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_BLACK**](stdfunc_8h.md#define-color_bright_black)  `(isatty(STDOUT\_FILENO) == true ? "\033[90m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_BLUE**](stdfunc_8h.md#define-color_bright_blue)  `(isatty(STDOUT\_FILENO) == true ? "\033[94m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_CYAN**](stdfunc_8h.md#define-color_bright_cyan)  `(isatty(STDOUT\_FILENO) == true ? "\033[96m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_GREEN**](stdfunc_8h.md#define-color_bright_green)  `(isatty(STDOUT\_FILENO) == true ? "\033[92m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_MAGENTA**](stdfunc_8h.md#define-color_bright_magenta)  `(isatty(STDOUT\_FILENO) == true ? "\033[95m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_RED**](stdfunc_8h.md#define-color_bright_red)  `(isatty(STDOUT\_FILENO) == true ? "\033[91m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_WHITE**](stdfunc_8h.md#define-color_bright_white)  `(isatty(STDOUT\_FILENO) == true ? "\033[97m" : "")`<br> |
| define  | [**COLOR\_BRIGHT\_YELLOW**](stdfunc_8h.md#define-color_bright_yellow)  `(isatty(STDOUT\_FILENO) == true ? "\033[93m" : "")`<br> |
| define  | [**COLOR\_CYAN**](stdfunc_8h.md#define-color_cyan)  `(isatty(STDOUT\_FILENO) == true ? "\033[36m" : "")`<br> |
| define  | [**COLOR\_DEFAULT**](stdfunc_8h.md#define-color_default)  `(isatty(STDOUT\_FILENO) == true ? "\033[0m" : "")`<br> |
| define  | [**COLOR\_ERROR**](stdfunc_8h.md#define-color_error)  `COLOR\_RED`<br> |
| define  | [**COLOR\_GREEN**](stdfunc_8h.md#define-color_green)  `(isatty(STDOUT\_FILENO) == true ? "\033[32m" : "")`<br> |
| define  | [**COLOR\_INFO**](stdfunc_8h.md#define-color_info)  `COLOR\_GREEN`<br> |
| define  | [**COLOR\_MAGENTA**](stdfunc_8h.md#define-color_magenta)  `(isatty(STDOUT\_FILENO) == true ? "\033[35m" : "")`<br> |
| define  | [**COLOR\_PURPLE**](stdfunc_8h.md#define-color_purple)  `(isatty(STDOUT\_FILENO) == true ? "\033[35m" : "")`<br> |
| define  | [**COLOR\_RED**](stdfunc_8h.md#define-color_red)  `(isatty(STDOUT\_FILENO) == true ? "\033[31m" : "")`<br> |
| define  | [**COLOR\_WARNING**](stdfunc_8h.md#define-color_warning)  `COLOR\_PURPLE`<br> |
| define  | [**COLOR\_WHITE**](stdfunc_8h.md#define-color_white)  `(isatty(STDOUT\_FILENO) == true ? "\033[37m" : "")`<br> |
| define  | [**COLOR\_YELLOW**](stdfunc_8h.md#define-color_yellow)  `(isatty(STDOUT\_FILENO) == true ? "\033[33m" : "")`<br> |
| define  | [**ERROR**](stdfunc_8h.md#define-error) (info) `{std::cout&lt;&lt;"--  ["&lt;&lt;COLOR\_RED&lt;&lt;"Error"&lt;&lt;COLOR\_DEFAULT&lt;&lt;"]: "&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl; exit(0);}`<br> |
| define  | [**ERROR\_COUT**](stdfunc_8h.md#define-error_cout)  `(isatty(STDOUT\_FILENO) == true ? "[\033[31mError: \033[0m] " : "[Error: ]")`<br> |
| define  | [**STATUS**](stdfunc_8h.md#define-status) (info) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_color**](stdfunc_8h.md#define-status_color) (info, color) `std::cout&lt;&lt;"--  "&lt;&lt;color&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_system\_time**](stdfunc_8h.md#define-status_system_time) (info, duration) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;", time: "&lt;&lt; double(duration) \* std::chrono::microseconds::period::num / std::chrono::microseconds::period::den&lt;&lt;" s"&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**STATUS\_time**](stdfunc_8h.md#define-status_time) (info, time\_taken) `std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_GREEN&lt;&lt;info&lt;&lt;", time: "&lt;&lt;(double)(time\_taken)/CLOCKS\_PER\_SEC&lt;&lt;" s"&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;`<br> |
| define  | [**WAIT**](stdfunc_8h.md#define-wait) (where) `{std::cout&lt;&lt;"Waiting "&lt;&lt;where&lt;&lt; ". Enter to continue..." &lt;&lt; std::endl; std::string dummy; std::getline(std::cin, dummy);}`<br> |
| define  | [**WARNING**](stdfunc_8h.md#define-warning) (info) `{std::cout&lt;&lt;"--  "&lt;&lt;COLOR\_YELLOW&lt;&lt;info&lt;&lt;COLOR\_DEFAULT&lt;&lt;std::endl;}`<br> |
| define  | [**WARN\_COUT**](stdfunc_8h.md#define-warn_cout)  `(isatty(STDOUT\_FILENO) == true ? "[\033[33mWarning: \033[0m] " : "[Warning: ]")`<br> |
| define  | [**xTHERMO\_VAR**](stdfunc_8h.md#define-xthermo_var)  <br> |

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
#define COLOR_BLUE `(isatty(STDOUT_FILENO) == true ? "\033[34m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_BLACK 

```C++
#define COLOR_BRIGHT_BLACK `(isatty(STDOUT_FILENO) == true ? "\033[90m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_BLUE 

```C++
#define COLOR_BRIGHT_BLUE `(isatty(STDOUT_FILENO) == true ? "\033[94m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_CYAN 

```C++
#define COLOR_BRIGHT_CYAN `(isatty(STDOUT_FILENO) == true ? "\033[96m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_GREEN 

```C++
#define COLOR_BRIGHT_GREEN `(isatty(STDOUT_FILENO) == true ? "\033[92m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_MAGENTA 

```C++
#define COLOR_BRIGHT_MAGENTA `(isatty(STDOUT_FILENO) == true ? "\033[95m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_RED 

```C++
#define COLOR_BRIGHT_RED `(isatty(STDOUT_FILENO) == true ? "\033[91m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_WHITE 

```C++
#define COLOR_BRIGHT_WHITE `(isatty(STDOUT_FILENO) == true ? "\033[97m" : "")`
```




<hr>



### define COLOR\_BRIGHT\_YELLOW 

```C++
#define COLOR_BRIGHT_YELLOW `(isatty(STDOUT_FILENO) == true ? "\033[93m" : "")`
```




<hr>



### define COLOR\_CYAN 

```C++
#define COLOR_CYAN `(isatty(STDOUT_FILENO) == true ? "\033[36m" : "")`
```




<hr>



### define COLOR\_DEFAULT 

```C++
#define COLOR_DEFAULT `(isatty(STDOUT_FILENO) == true ? "\033[0m" : "")`
```




<hr>



### define COLOR\_ERROR 

```C++
#define COLOR_ERROR `COLOR_RED`
```




<hr>



### define COLOR\_GREEN 

```C++
#define COLOR_GREEN `(isatty(STDOUT_FILENO) == true ? "\033[32m" : "")`
```




<hr>



### define COLOR\_INFO 

```C++
#define COLOR_INFO `COLOR_GREEN`
```




<hr>



### define COLOR\_MAGENTA 

```C++
#define COLOR_MAGENTA `(isatty(STDOUT_FILENO) == true ? "\033[35m" : "")`
```




<hr>



### define COLOR\_PURPLE 

```C++
#define COLOR_PURPLE `(isatty(STDOUT_FILENO) == true ? "\033[35m" : "")`
```




<hr>



### define COLOR\_RED 

```C++
#define COLOR_RED `(isatty(STDOUT_FILENO) == true ? "\033[31m" : "")`
```




<hr>



### define COLOR\_WARNING 

```C++
#define COLOR_WARNING `COLOR_PURPLE`
```




<hr>



### define COLOR\_WHITE 

```C++
#define COLOR_WHITE `(isatty(STDOUT_FILENO) == true ? "\033[37m" : "")`
```




<hr>



### define COLOR\_YELLOW 

```C++
#define COLOR_YELLOW `(isatty(STDOUT_FILENO) == true ? "\033[33m" : "")`
```




<hr>



### define ERROR 

```C++
#define ERROR (
    info
) `{std::cout<<"--  ["<<COLOR_RED<<"Error"<<COLOR_DEFAULT<<"]: "<<info<<COLOR_DEFAULT<<std::endl; exit(0);}`
```




<hr>



### define ERROR\_COUT 

```C++
#define ERROR_COUT `(isatty(STDOUT_FILENO) == true ? "[\033[31mError: \033[0m] " : "[Error: ]")`
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
#define WARN_COUT `(isatty(STDOUT_FILENO) == true ? "[\033[33mWarning: \033[0m] " : "[Warning: ]")`
```




<hr>



### define xTHERMO\_VAR 

```C++
#define xTHERMO_VAR 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/stdfunc.h`

