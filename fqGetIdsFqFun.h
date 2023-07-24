/*#########################################################
# Name: fqGetIdsFqFun
# Use:
#  - Holds functions used to process the fastq header,
#    move to the next fastq entry, and print out fastq
#    entries.
# Libraries:
#  - "vectorWrap.h"
#  - "fqGetIdsStructs.h"
# C Standard Libraries
#  o <limits.h>
#  o <stdint.h>
#  o <stdio.h>
#  o <stdlib.h>
#  o <immintrin.h>
#########################################################*/

#ifndef FQGREPFQFUN_H
#define FQGREPFQFUN_H

#include "vectorWrap.h"  // Using vectors
#include "fqGetIdsStructs.h" /*<stdlib.h>, <stdio.h>, <stdint.h>*/
     /*Has structure functions & hex look up tables*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
' o fun-01 parseFastqHeader:
'   - Reads in and converts read id in buffer to al number
' o fun-02/fun-04 moveToNextFastqEntry:
Usesing'   - Moves the the next fastq entry and prints to the
'     output file if an output file was provided.
'   - fun-02 is the vector version of moveToNextFastqEntry
'   - fun-04 has the same header and is the scaler version
'     moveToNextFastqEntry
' o fun-03 repositionBuff:
'   - Reposition posUL to point to the element you want
'     to move past in the buffer (true value in result
'     vect).
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------/
| Output:
|  - Modifies:
|    o startNameCStr to point to the start of the read name
|    o posUL to point to the end of the read id
|    o buffCStr to hold more input if needed
|    o LenIdInt to hold the length of the read id if more
|      input was raead in
|  - Returns:
|    o 0 if nothing went wrong
|    o 1 If the end of the file
|    o 4 If ran out of file
\--------------------------------------------------------*/
uint8_t parseFastqHeader(
    char *buffCStr,          // Buffer with reads
    unsigned long *posUL,    // Point working at in buffer
    //unsigned long *startUL,
       // Points to start of read entry. This is changed
       // when new input is grabed from fastqFile
    unsigned long *lenInUL,  //Number of char in buffer
    unsigned int buffSizeUI, // Size of the buffer
    struct bigNum *idBigNum, // Will hold big number
    FILE *fqFile             // Fastq file to get data from
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: parseFastqHeader
   '  - Reads in and converts read id in buffer to a
   '    number.
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Move to start of read id
   '  o fun-01 sec-03:
   '    - Convert read id to bigNum read id (long array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Prints:
|    o Reads to outFILE if outFILE != 0
|  - Modifies:
|    o buffCStr to have the next buffer if empty
|    o incurments pointInBufferCStr to start of next read
|  - Returns:
|   o 0 if nothing went wrong
|   o 1 If the end of the file
|   o 4 If ran out of file
\--------------------------------------------------------*/
uint8_t moveToNextFastqEntry(
    char *buffCStr,       // buffer with input to scan
    unsigned long *posUL, // Position at in buffer
    unsigned int lenBuffI,// Size of buffer to work on
    unsigned long *numReadCharUL,//Number chars in buffCStr
    FILE *fqFILE,         // Fastq file to get input from
    FILE *outFILE         // FILE to output reads to
      // This function will skip printing output if
      // outFILE == 0
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: or Fun-04 TOC: moveToNextFastqEntry
   '  - Moves the the next fastq entry and prints to the
   '    output file if an output file was provided.
   '  - This is the vector version of moveToNextFastqEntry
   '  o fun-02 sec-01: or fun-04 sec-01:
   '    -  Variable declerations
   '  o fun-02 sec-02: or fun-04 sec-02:
   '    - Move past header
   '  o fun-02 sec-03: or fun-04 sec-03:
   '    - Find number newlines in sequence & move to spacer
   '  o fun-02 sec-04: or fun-04 sec-04:
   '    - Move past spacer
   '  o fun-02 Sec-05: or fun-04 sec-05:
   '    - Move past Q-score entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#endif
