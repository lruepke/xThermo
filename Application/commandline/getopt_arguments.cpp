#include "getopt_arguments.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <string.h>
#endif

#include <stdlib.h>

char* optarg = NULL;
int optind = 2;

// An argument that can be the value of an option: anything that does not start with '-', or a number such as -0.1
static int is_value(const char *arg)
{
    char *end = NULL;
    if (arg[0] != '-') return 1;
    strtod(arg, &end);
    return end != arg && *end == 0;
}

// Options are single letters. A value follows the letter directly (-P300E5) or as the next argument (-P 300E5).
// Returns the letter, '?' for an unknown option, ':' for an option whose value is missing, and -1 at the end of
// the options (no more arguments, or an argument that is not an option).
int getopt_arguments(int argc, char *const argv[], const char *optstring)
{
    if (optind >= argc || argv[optind][0] != '-' || argv[optind][1] == 0)
    {
        return -1;
    }
    int opt = argv[optind][1];
    const char *p = strchr(optstring, opt);
    if (p == NULL || opt == ':')
    {
        optind++;
        return '?';
    }
    if (p[1] == ':') //option with a value
    {
        if (argv[optind][2] != 0)
        {
            optarg = argv[optind] + 2;
            optind++;
        } else if (optind + 1 < argc && is_value(argv[optind + 1]))
        {
            optarg = argv[optind + 1];
            optind += 2;
        } else
        {
            optind++;
            return ':';
        }
    } else //option without a value: nothing may follow the letter, except for -h and -v (e.g. -help)
    {
        if (argv[optind][2] != 0 && opt != 'h' && opt != 'v')
        {
            optind++;
            return '?';
        }
        optind++;
    }

    return opt;
}
