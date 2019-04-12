/*
 * This is from data.table devel version 1.10.5
 * Slightly modified, or heavily, however you wish to rank that.
 * Huge thanks to Matt Dowle for his help and willingness to let me play
 * 
 */

#ifndef FWRITE
#define FWRITE


#define STRICT_R_HEADERS
#define STOP     error
#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdint.h>

typedef void (*writer_fun_t)(void *, int64_t, char **);

void writeInt32();
void writeFloat64();
void writeString();

void write_chars(const char *source, char **dest);

SEXP fwriteMain(SEXP MAT, SEXP filename_Arg, SEXP sep_Arg,
                SEXP eol_Arg, SEXP dec_Arg, SEXP buffMB_Arg);


#endif