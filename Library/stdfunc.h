#ifndef STDFUNC_H
#define STDFUNC_H
#include <cstdio>
#include <iostream>
#include<map>
#include <vector>
#include <sys/stat.h>
#include <cstring>
#include <sstream>
#include <iomanip> //to use std::setprecision
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#ifndef xTHERMO_VAR

#ifdef WIN32
#  ifdef xTHERMO_DLL
#    ifdef xThermal_DLL_EXPORT
#      define xTHERMO_VAR __declspec(dllexport)
#    else
#      define xTHERMO_VAR __declspec(dllimport)
#    endif
#  else
#    define xTHERMO_VAR
#  endif
#else
#  define xTHERMO_VAR
#endif

#endif

#ifdef _WIN32
    // define color, this seems only work on MacOS and linux, doesn't work on windows
    #include "windows.h"
    #define BLACK			0
    #define BLUE			1
    #define GREEN			2
    #define CYAN			3
    #define RED				4
    #define MAGENTA			5
    #define BROWN			6
    #define LIGHTGRAY		7
    #define DARKGRAY		8
    #define LIGHTBLUE		9
    #define LIGHTGREEN		10
    #define LIGHTCYAN		11
    #define LIGHTRED		12
    #define LIGHTMAGENTA	13
    #define YELLOW			14
    #define WHITE			15
    #define BRIGHT_BLACK	0
    #define BRIGHT_RED	    4
    #define BRIGHT_GREEN	2
    #define BRIGHT_YELLOW	14
    #define BRIGHT_BLUE		1
    #define BRIGHT_MAGENTA	5
    #define BRIGHT_CYAN		3
    #define BRIGHT_WHITE	15

    static HANDLE   m_hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
    static WORD     m_currentConsoleAttr;
    static CONSOLE_SCREEN_BUFFER_INFO csbi;

    #define ERROR_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<"Error: "<<COLOR_DEFAULT
    #define WARN_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<"Warning: "<<COLOR_DEFAULT
    #define COLOR_PURPLE ""
    #define COLOR_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<""
    #define COLOR_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (GREEN & 0x0F) );cout<<""
    #define COLOR_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<""
    #define COLOR_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BLUE & 0x0F) );cout<<""
    #define COLOR_MAGENTA "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (MAGENTA & 0x0F) );cout<<""
    #define COLOR_CYAN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (CYAN & 0x0F) );cout<<""
    #define COLOR_WHITE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (WHITE & 0x0F) );cout<<""
    #define COLOR_BRIGHT_BLACK "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_BLACK & 0x0F) );cout<<""
    #define COLOR_BRIGHT_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_RED & 0x0F) );cout<<""
    #define COLOR_BRIGHT_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_GREEN & 0x0F) );cout<<""
    #define COLOR_BRIGHT_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_YELLOW & 0x0F) );cout<<""
    #define COLOR_BRIGHT_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_BLUE & 0x0F) );cout<<""
    #define COLOR_BRIGHT_MAGENTA "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_MAGENTA & 0x0F) );cout<<""
    #define COLOR_BRIGHT_CYAN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_CYAN & 0x0F) );cout<<""
    #define COLOR_BRIGHT_WHITE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BRIGHT_WHITE & 0x0F) );cout<<""
    #define COLOR_DEFAULT "";SetConsoleTextAttribute(m_hConsole, m_currentConsoleAttr );cout<<""
#else
    #include <unistd.h> // use isatty(STDOUT_FILENO) check if console is used or redirect to file, then decide use colored text or not
    // define color, this seems only work on MacOS and linux, doesn't work on windows: see https://en.wikipedia.org/wiki/ANSI_escape_code
    #define ERROR_COUT (isatty(STDOUT_FILENO) == true ? "[\033[31mError: \033[0m] " : "[Error: ]")
    #define WARN_COUT (isatty(STDOUT_FILENO) == true ? "[\033[33mWarning: \033[0m] " : "[Warning: ]")
    #define COLOR_PURPLE (isatty(STDOUT_FILENO) == true ? "\033[35m" : "")
    #define COLOR_RED (isatty(STDOUT_FILENO) == true ? "\033[31m" : "")
    #define COLOR_GREEN (isatty(STDOUT_FILENO) == true ? "\033[32m" : "")
    #define COLOR_YELLOW (isatty(STDOUT_FILENO) == true ? "\033[33m" : "")
    #define COLOR_BLUE (isatty(STDOUT_FILENO) == true ? "\033[34m" : "")
    #define COLOR_MAGENTA (isatty(STDOUT_FILENO) == true ? "\033[35m" : "")
    #define COLOR_CYAN (isatty(STDOUT_FILENO) == true ? "\033[36m" : "")
    #define COLOR_WHITE (isatty(STDOUT_FILENO) == true ? "\033[37m" : "")
    #define COLOR_BRIGHT_BLACK (isatty(STDOUT_FILENO) == true ? "\033[90m" : "")
    #define COLOR_BRIGHT_RED (isatty(STDOUT_FILENO) == true ? "\033[91m" : "")
    #define COLOR_BRIGHT_GREEN (isatty(STDOUT_FILENO) == true ? "\033[92m" : "")
    #define COLOR_BRIGHT_YELLOW (isatty(STDOUT_FILENO) == true ? "\033[93m" : "")
    #define COLOR_BRIGHT_BLUE (isatty(STDOUT_FILENO) == true ? "\033[94m" : "")
    #define COLOR_BRIGHT_MAGENTA (isatty(STDOUT_FILENO) == true ? "\033[95m" : "")
    #define COLOR_BRIGHT_CYAN (isatty(STDOUT_FILENO) == true ? "\033[96m" : "")
    #define COLOR_BRIGHT_WHITE (isatty(STDOUT_FILENO) == true ? "\033[97m" : "")
    #define COLOR_DEFAULT (isatty(STDOUT_FILENO) == true ? "\033[0m" : "")
#endif

namespace xThermal{

    #define STATUS(info) std::cout<<"--  "<<COLOR_GREEN<<info<<COLOR_DEFAULT<<std::endl;
    #define STATUS_color(info, color) std::cout<<"--  "<<color<<info<<COLOR_DEFAULT<<std::endl;
    #define STATUS_time(info, time_taken) std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<<(double)(time_taken)/CLOCKS_PER_SEC<<" s"<<COLOR_DEFAULT<<std::endl;
    #define STATUS_system_time(info, duration) std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<< double(duration) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den<<" s"<<COLOR_DEFAULT<<std::endl;
    #define ERROR(info) {std::cout<<"--  ["<<COLOR_RED<<"Error"<<COLOR_DEFAULT<<"]: "<<info<<COLOR_DEFAULT<<std::endl; exit(0);}
    #define WARNING(info) {std::cout<<"--  "<<COLOR_YELLOW<<info<<COLOR_DEFAULT<<std::endl;}
    #define ASSERT(expression, info) {if(!(expression))std::cout<<"--  "<<COLOR_RED<<info<<COLOR_DEFAULT<<std::endl; assert(expression);}
    #define WAIT(where) {std::cout<<"Waiting "<<where<< ". Enter to continue..." << std::endl; std::string dummy; std::getline(std::cin, dummy);}
    // some special definition
    #define COLOR_INFO COLOR_GREEN
    #define COLOR_ERROR COLOR_RED
    #define COLOR_WARNING COLOR_PURPLE

    //some utility functions
    std::vector<std::string> string_split (const std::string& s, const std::string& delimiter) ;
    std::string extname_file(const std::string& filepath);
    std::string filename_without_ext(const std::string& filepath);

}
#endif