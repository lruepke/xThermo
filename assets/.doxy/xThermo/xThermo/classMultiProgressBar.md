

# Class MultiProgressBar



[**ClassList**](annotated.md) **>** [**MultiProgressBar**](classMultiProgressBar.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MultiProgressBar**](#function-multiprogressbar-12) (std::vector&lt; double &gt; xmin, std::vector&lt; double &gt; xmax, std::vector&lt; std::string &gt; title) <br> |
|   | [**MultiProgressBar**](#function-multiprogressbar-22) (double total, int color=0) <br> |
|  void | [**Update**](#function-update-12) (double current\_pos=-1) <br> |
|  void | [**Update**](#function-update-22) (std::vector&lt; double &gt; current\_pos) <br> |
|   | [**~MultiProgressBar**](#function-multiprogressbar) () <br> |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  char | [**m\_bar\_char\_left**](#variable-m_bar_char_left)  <br> |
|  char | [**m\_bar\_char\_right**](#variable-m_bar_char_right)  <br> |
|  int | [**m\_bar\_number**](#variable-m_bar_number)  <br> |
|  std::vector&lt; std::string &gt; | [**m\_bar\_str**](#variable-m_bar_str)  <br> |
|  std::vector&lt; double &gt; | [**m\_current\_index**](#variable-m_current_index)  <br> |
|  std::vector&lt; double &gt; | [**m\_left**](#variable-m_left)  <br> |
|  int | [**m\_length\_bar**](#variable-m_length_bar)  <br> |
|  int | [**m\_maxLength\_title**](#variable-m_maxlength_title)  <br> |
|  std::vector&lt; double &gt; | [**m\_percent**](#variable-m_percent)  <br> |
|  std::vector&lt; double &gt; | [**m\_right**](#variable-m_right)  <br> |
|  std::vector&lt; std::string &gt; | [**m\_title**](#variable-m_title)  <br> |
|  std::vector&lt; double &gt; | [**m\_total**](#variable-m_total)  <br> |




















## Public Functions Documentation




### function MultiProgressBar [1/2]

```C++
MultiProgressBar::MultiProgressBar (
    std::vector< double > xmin,
    std::vector< double > xmax,
    std::vector< std::string > title
) 
```




<hr>



### function MultiProgressBar [2/2]

```C++
MultiProgressBar::MultiProgressBar (
    double total,
    int color=0
) 
```




<hr>



### function Update [1/2]

```C++
void MultiProgressBar::Update (
    double current_pos=-1
) 
```




<hr>



### function Update [2/2]

```C++
void MultiProgressBar::Update (
    std::vector< double > current_pos
) 
```




<hr>



### function ~MultiProgressBar 

```C++
MultiProgressBar::~MultiProgressBar () 
```




<hr>
## Protected Attributes Documentation




### variable m\_bar\_char\_left 

```C++
char MultiProgressBar::m_bar_char_left;
```




<hr>



### variable m\_bar\_char\_right 

```C++
char MultiProgressBar::m_bar_char_right;
```




<hr>



### variable m\_bar\_number 

```C++
int MultiProgressBar::m_bar_number;
```




<hr>



### variable m\_bar\_str 

```C++
std::vector<std::string> MultiProgressBar::m_bar_str;
```




<hr>



### variable m\_current\_index 

```C++
std::vector<double> MultiProgressBar::m_current_index;
```




<hr>



### variable m\_left 

```C++
std::vector<double> MultiProgressBar::m_left;
```




<hr>



### variable m\_length\_bar 

```C++
int MultiProgressBar::m_length_bar;
```




<hr>



### variable m\_maxLength\_title 

```C++
int MultiProgressBar::m_maxLength_title;
```




<hr>



### variable m\_percent 

```C++
std::vector<double> MultiProgressBar::m_percent;
```




<hr>



### variable m\_right 

```C++
std::vector<double> MultiProgressBar::m_right;
```




<hr>



### variable m\_title 

```C++
std::vector<std::string> MultiProgressBar::m_title;
```




<hr>



### variable m\_total 

```C++
std::vector<double> MultiProgressBar::m_total;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/MultiProgressBar.h`

