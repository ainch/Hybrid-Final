#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "def.h"

#ifndef __XYPIC_H__
#define __XYPIC_H__
// C Function declaration
void display_title();
void display_finish();
void InputRead(int argc, char *argv[]);
void start();
void DumpRead(int argc, char *argv[]);
void MAKE_Value();
void MAKE_TECPLOT();
#endif