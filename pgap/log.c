/******************************************************************************
 * log.h: logging support                                                     *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef LOG_C
#define LOG_C

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include "log.h"

char glogPrefix[256] = {'g','a'};
char glogDir[256] = {'.'};
IDTYPE glogId = 0;

int debug = 0;

void glog_init(char *d, char *p, IDTYPE i)
{
    sprintf(glogDir, "%s", d);
    sprintf(glogPrefix, "%s", p);
    glogId = i;
}
// logId is usually process/thread id
void glog(char *format, ...)
{
    static char logBuf[LOG_BUF_SIZE];
    static FILE * logF = NULL; // file desc
    static time_t prev_t = 0; // last time log entry was recorded
    static char * logPtr = NULL; // pointer to current empty memory
    static char logtext[LOG_ENTRY_BUF_SIZE]; // log entry buffer
    static char fn[256]; // log file name
    
    // first time: set log path
    if (logPtr == NULL) {
        char hostname[128];
        if (gethostname(hostname, 128) == -1) {
            strcpy(hostname, "localhost");
        }
        //sprintf(fn, "%s/%s_%s_%d.log", logDir, logPrefix, hostname, logId);
        //TODO: to be portable to Windows: use '\'
        if (strchr(glogPrefix, '/'))
            sprintf(fn, "%s_%ld.log", glogPrefix, (long int)glogId);
        else
            sprintf(fn, "%s/%s_%ld.log", glogDir, glogPrefix, (long int)glogId);
        logPtr = logBuf;
        // init file
        if ((logF = fopen(fn, "a")) == NULL) {
            perror("Error opening log file");
            logF = stdout;
            return;
        }
    }
    // write log entry
    time_t t = time(NULL);
    char time_str[64];
    //strftime(time_str, 64, "%m/%d/%Y %T", gmtime(&t));
    strftime(time_str, 64, "%T", localtime(&t));
    va_list argp;
    va_start(argp, format);
    vsnprintf(logtext, LOG_ENTRY_BUF_SIZE, format, argp);
    va_end(argp);
    int len1, len2;
    len1 = strlen(time_str); len2 = strlen(logtext);
    if ((strlen(logBuf)+strlen(time_str)+strlen(logtext) > LOG_BUF_SIZE) \
       || (t - prev_t > 20) || (glogId == -1)) {
        // write buffer to file
        if (fprintf(logF, "%s", logBuf) < 0) {
            fprintf(stderr, "Error writing to log file %s\n", fn);
            fflush(stderr);
            if (logF != stdout)
                fclose(logF);
            logPtr = NULL;
        } else {
            fflush(logF);
            // reset pointer
            logPtr = logBuf; *logPtr='\0'; 
            if (glogId == -1) {
                fprintf(logF, "%s %s", time_str, logtext); 
                fclose(logF); // close log file
                return;
            }
        }
        prev_t = t;
    }
    sprintf(logPtr, "%s %s", time_str, logtext);
    //output to stderr for job scheduler execution
    //so that when walltime is up, we don't lose info
    fprintf(stderr, "@@%s %s", time_str, logtext);
    logPtr +=  len1 + 1 + len2 + 1;
    prev_t = t;
}

#endif
