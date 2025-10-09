#ifndef GETOPT_H__
#define GETOPT_H__

#ifdef __cplusplus
extern "C" {
#endif

extern char *optarg;
extern int optind;

int getopt_arguments(int argc, char *const argv[], const char *optstring);

#ifdef __cplusplus
}
#endif

#endif