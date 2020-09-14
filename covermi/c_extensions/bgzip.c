/* bgzip.c -- Block compression/decompression utility.

   Copyright (C) 2008, 2009 Broad Institute / Massachusetts Institute of Technology
   Copyright (C) 2010, 2013-2018 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notices and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/stat.h>
#include "bgzf.h"


static const int WINDOW_SIZE = 64 * 1024;


int main(int argc, char **argv) {
    int c;
    BGZF *fp;
    void *buffer;
    long start, end;

    fp = bgzf_open("fn", "r");
    if (fp == NULL) {
        fprintf(stderr, "[bgzip] Could not open file: %s\n");
        return 1;
        }

    if (!fp->is_compressed) {
        fprintf(stderr, "[bgzip] Expected compressed file -- ignored\n");
        return 1;
        }

    buffer = malloc(WINDOW_SIZE);
    
    if ( start>0 ) {
        if ( bgzf_index_load(fp, argv[optind], ".gzi") < 0 ) {
            error("Could not load index: %s.gzi\n", argv[optind]);
            }
        if ( bgzf_useek(fp, start, SEEK_SET) < 0 ) {
            error("Could not seek to %d-th (uncompressd) byte\n", start);
            }
        }

    while (1) {
        if (end < 0) {
            c = bgzf_read(fp, buffer, WINDOW_SIZE);
            }
        else { 
            c = bgzf_read(fp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
            }
        if (c == 0) {
            break;
            }
        if (c < 0) {
            error("Error %d in block starting at offset %" PRId64 "(%" PRIX64 ")\n", fp->errcode, fp->block_address, fp->block_address);
            }
        start += c;
        if (end >= 0 && start >= end) {
            break;
            }
        }
    free(buffer);
    if (bgzf_close(fp) < 0) {
        error("Close failed: Error %d\n",fp->errcode);
        }
    return 0;
    }
