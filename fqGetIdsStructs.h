/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepStructs TOC:
#    struct-1: readInfo: Holds the read name
#    struct-2: readNodeStack: Makes a stack (filo) of read nodes
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef FQGREPSTRUCTS_H
#define FQGREPSTRUCTS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h> /*for intx_t & uintx_t variables*/

/*Look up table to use in converting char to hext*/
extern char hexTblCharAry[];

#define defBitsPerChar 5
#define defBitsInL (sizeof(long) << 3) - 1 /*-1 to account for - flag*/
#define defBitsInUL (sizeof(unsigned long) << 3)

#define defBitsPerLimb (sizeof(long) << 3) - 1 /*-1 for - flag*/
#define defMaxDigPerLimb (((sizeof(long) << 3)-1)/defBitsPerChar)

/*Settings for memory efficent or speed mode*/
/*defOSBit is to mark if 32bit system or not*/


#if INTPTR_MAX >= INT64_MAX
  #define defOSBit 64
#elif INTPTR_MAX >= INT32_MAX
  #define defOSBit 32
#elif INTPTR_MAX >= INT16_MAX
  #define defOSBit 16
#else
  #error System does not support 16 bit integers
#endif

/*---------------------------------------------------------------------\
| Struct-3: bigNum
| Use: Stores hex components of string as a big number (array of longs)
\---------------------------------------------------------------------*/
typedef struct bigNum
{ /*bigNum*/
    /*Memory (MEM) or speed version (default)*/
    long totalL;                  /*All limbs added together*/
    long *bigNumAryIOrL;

    unsigned char lenUsedElmChar; /*Number of elements used long array*/
    unsigned char lenAllElmChar;   /*Size of the long array*/
}bigNum;

/*---------------------------------------------------------------------\
| Struct-1: readInfo
| Use: Stores the hex id of a read as a big number
\---------------------------------------------------------------------*/
typedef struct readInfo
{ /*readInfo structer*/
    int8_t balanceChar;     /*Tells if the node is balanced*/
    struct bigNum *idBigNum; /*Holds read id as unique big number*/
    struct readInfo *leftChild; 
    struct readInfo *rightChild;
}readInfo; /*readInfo structure*/

/*---------------------------------------------------------------------\
| Output: Modifies: readInfoStruct to have default values (all 0's)
\---------------------------------------------------------------------*/
struct readInfo * makeReadInfoStruct(
    char *readIdCStr,         /*c-string with read name to copy*/
    const unsigned int *lenCStrUInt /*Length of cString to convert*/
); /*Allocates memomory & makes a readInfo structer (set to 0)*/

void freeReadInfoStruct(
    struct readInfo **readInfoStruct /*struct to free*/
); /*frees a readInfo structer*/

/*---------------------------------------------------------------------\
| Struct-2: readNodeStack
| Use: Used to make a stack of readInfo structers, used in transvering my
|      tree of reads
| Note:
|   o Users makes array of nodes they pass around
|   o A well balanced AVL tree coverts 10^18 nodes at a depth of 73
\---------------------------------------------------------------------*/
typedef struct readNodeStack
{ /*readNodeStack*/
    struct readInfo *readNode;                /*Node in stack*/
}readNodeStack; /*readNodeStack*/

/*---------------------------------------------------------------------\
| Output:
|    Modifes readStack to piont to last element
\---------------------------------------------------------------------*/
void pushReadNodeStack(
    struct readNodeStack **readStackAry, /*Array of read info nodes*/
    struct readInfo *readNode        /*readInfo to assing to next node*/
); /*pushes a readNodeStack structer onto a readNodeStack stack*/

/*---------------------------------------------------------------------\
| Output: Modifies: readNodeStack to point to next node in stack
\---------------------------------------------------------------------*/
void popReadNodeStack(
    struct readNodeStack **readStackAry /*readInfo Array (stack) to pop*/
); /*sets readNodeStack to next readInfo node in stack*/

/*---------------------------------------------------------------------\
| Output:
|    returns: bigNum structer with the converted big number
|    returns: 0 if memory alloaction failed
\---------------------------------------------------------------------*/
struct bigNum * makeBigNumStruct(
    char *cStrToCnvt,/*C-string to convert hex elements to big number*/
    const unsigned int *lenCStrUInt /*Length of cString to convert*/
); /*Convert hex characters in c-string to bigNum struct*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies: idBigNum to hold the converted big number.
|        - Sets idBigNum->lenUsedULng to 0 if memory reallocation failed
\---------------------------------------------------------------------*/
void strToBackwardsBigNum(
    struct bigNum *idBigNum,  /*Holds the output big number*/
    char *cStrToCnvt,         /*C-string to convert to large number*/
    const unsigned int *lenCStrUInt /*Length of cString to convert*/
); /*Flips c-string & converts to big number*/

/*---------------------------------------------------------------------\
| Output:
|    frees: idBigNum & sets pointer to 0
\---------------------------------------------------------------------*/
void freeBigNumStruct(
    struct bigNum **idBigNum  /*Address to bigNum structer to free*/
); /*Frees a bigNum struct from memory*/

/*---------------------------------------------------------------------\
| Output:
|    Returns:
|      - 0 if both equal,
|      - < 0 if query is smaller than subject
|      - > 0 if qeury is bigger than subject
\---------------------------------------------------------------------*/
long cmpBigNums(
    struct bigNum *bigNumOne,
    struct bigNum *bigNumTwo
); /*Compares bigNumOne to bigNumTwo to see if equal, >, or <*/

/*---------------------------------------------------------------------\
| Output:
|    Modifies:
|      - EndNameCStr to point to '\n', ' ', & '\t' at end of read name
|      - LenIdULng to hold the length of the read id
|      - LenInputInt to hold length of input buffer
|    Returns:
|      - 0 if fails or end of file (lenIdULng < buffSizeInt)
|      - pointer to struct with bigNum struct having converted read id
\---------------------------------------------------------------------*/
struct readInfo * cnvtIdToBigNum(
    char *bufferCStr, /*buffer to hold fread input (can have data)*/
    uint32_t buffSizeInt,         /*Size of buffer to work on*/
    char **endNameCStr, /*Points to start of id, will point to end*/
    uint64_t *lenInputULng,        /*Length of input from fread*/
    unsigned char *lenBigNumChar, /*Holds size to make bigNumber*/
    FILE *idFILE          /*Fastq file to get data from*/
); /*Converts read id to bigNum read id, will grab new file input*/

#endif
