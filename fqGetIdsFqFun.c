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

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
' o fun-01 parseFastqHeader:
'   - Reads in and converts read id in buffer to al number
' o fun-02 moveToNextFastqEntry:
'   - Moves the the next fastq entry and prints to the
'     output file if an output file was provided.
'   - This is the vector version of moveToNextFastqEntry
' o fun-04 moveToNextFastqEntry:
'   - Moves the the next fastq entry and prints to the
'     output file if an output file was provided.
'   - This is the no vector version of moveToNextFastqEntry
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "fqGetIdsFqFun.h"

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
    FILE *fqFILE             // Fastq file to get data from
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: parseFastqHeader
   '  - Reads in and converts read id in buffer to a
   '    number.
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Check if need to get more file input
   '  o fun-01 sec-03:
   '    - Convert read id to bigNum read id (long array)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t charBit = 0;

    /*For speed versus memroy (MEM) compile settings*/
    long *elmILPtr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-02:
    ^  - Check if need to get more file input
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    idBigNum->lenUsedElmChar = 0;
    idBigNum->totalL = 0;

    if(*posUL >= *lenInUL - 129) 
    { // If I need to get new input
        /*move though fastq function will catch EOF*/
	if(*lenInUL < buffSizeUI) return 4;

        fseek(fqFILE,(*lenInUL - *posUL)*-1,SEEK_CUR);
        *lenInUL =
          fread(buffCStr,sizeof(char),buffSizeUI,fqFILE);

        *(buffCStr + *lenInUL) = '\0';
        *posUL = 0;
    } // If I need to get new input
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-01 Sec-03:
    ^  - Convert read id to bigNum read id (long array)
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    ++(*posUL); /*Get off the @*/

     do { /*While still on the read name part of header*/
        /*Graph unsigned long element working on*/
        elmILPtr =
            idBigNum->bigNumAryIOrL
          + idBigNum->lenUsedElmChar;

        /*Track number used longs*/
        (idBigNum->lenUsedElmChar)++;
        *elmILPtr = 0;
        charBit = 0;

        while(charBit < defMaxDigPerLimb)
        { // while have room in the current long
            if(
              hexTblCharAry[
                (unsigned char) *(buffCStr + *posUL)
              ] & 64
            ) break; // finshed converting the hex string

            if(
              !(hexTblCharAry[
                 (unsigned char) *(buffCStr + *posUL)
               ] & 32)
            ) { /*If is a hex character*/
                *elmILPtr = *elmILPtr << defBitsPerChar;

                *elmILPtr +=
                  hexTblCharAry[
                    (unsigned char) *(buffCStr + *posUL)
                ]; // Convert the character to a number

                ++charBit;
            } /*If is a hex character*/

            ++(*posUL);
        } // while have room in the current long

        /*Only used in speed setting*/
        #if defOSBit == 64
          idBigNum->totalL +=
            *((int*) elmILPtr) + *(((int*) elmILPtr) + 1);
        #elif defOSBit == 32
          idBigNum->totalL +=
            *((short*) elmILPtr)+ *(((short*) elmILPtr)+1);
        #elif defOSBit == 16
          idBigNum->totalL +=
            *((char*) elmILPtr) + *(((char*) elmILPtr) +1);
        #else
          #error System does not support 16 bit integers
        #endif
    } while(
      !(hexTblCharAry[
         (unsigned char) *(buffCStr + *posUL)
      ] & 64)
    ); // While still on the read name part of header

     return 0; /*Copied name sucessfully*/
} /*parseFastqHeader*/

#if defined AVX512 || AVX2 || SSE

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: moveToNextFastqEntry
   '  - Moves the the next fastq entry and prints to the
   '    output file if an output file was provided.
   '  - This is the vector version of moveToNextFastqEntry
   '  o fun-02 sec-01:
   '    -  Variable declerations
   '  o fun-02 sec-02:
   '    - Move past header
   '  o fun-02 sec-03:
   '    - Find number newlines in sequence & move to spacer
   '  o fun-02 sec-04:
   '    - Move past spacer
   '  o fun-02 Sec-05:
   '    - Move past Q-score entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-02 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long startPosUL = *posUL;
    long numSeqLinesL = 0;

    vectI8 lineVect;
    vectI8 newLineVect; // For new line comparisons
    vectI8 plusVect;    // marker between seq and q-score

    mmaskI8 resMaskI8;
    mmaskI8 tmpMaskI8;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-02 Sec-02:
    ^  - Move past header
    ^  o fun-02 sec-02 sub-01:
    ^    - Initial set up / ensure have data in buffer
    ^  o fun-02 sec-02 sub-02:
    ^    - Move past the header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-02 Sec-02 Sub-01:
    ^  - Initial set up / ensure have data in buffer
    \*****************************************************/

    if(*numReadCharUL - *posUL < vectorBytes)
    { // IF I need to have an offset
      // Check if have more file to read in
      if(*numReadCharUL < lenBuffI) return 1;

      // Get back to char after '\n'
      fseek(
        fqFILE,
        -1 * (*numReadCharUL - *posUL),
        SEEK_CUR
      );

      *posUL = 0;
      startPosUL = 0;

      *numReadCharUL =
        fread(buffCStr, sizeof(char), lenBuffI, fqFILE);
    } // IF I need to have an offset

    mmLoadUI8(lineVect, (buffCStr + *posUL));
    mmSet1I8(newLineVect, '\n');
    mmSet1I8(plusVect, '+');

    mmCmpEqI8(resMaskI8, lineVect, newLineVect);

    /*****************************************************\
    * Fun-02 Sec-02 Sub-02:
    *  - Move past the header
    \*****************************************************/

    while(checkVectMaskTrueI8(resMaskI8) == 0)
    { /*While on header, move to sequence line*/

        if(*posUL + (vectorBytes << 1) >= *numReadCharUL)
        { /*If ran out of buffer */
            if(*numReadCharUL < lenBuffI)
                return 4;      // incomplete fastq entry

            if(outFILE != 0)
            { // If I am writing output to a file
              *(buffCStr + *posUL + vectorBytes) = '\0';
              fprintf(outFILE,"%s",buffCStr+startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *numReadCharUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file
            
            if(*posUL + vectorBytes < *numReadCharUL)
              /*Move to charaters that were not loaded*/
              fseek(
                fqFILE,
                -1 * (*numReadCharUL-*posUL-vectorBytes),
                SEEK_CUR
              );

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            // Keep as a c-string
            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            goto headerVectLoad;
        } /*If ran out of buffer*/

        *posUL += vectorBytes;

        headerVectLoad:

        mmLoadUI8(lineVect, (buffCStr + *posUL));
        mmCmpEqI8(resMaskI8, lineVect, newLineVect);
    } /*While on header, move to sequence line*/

    /*****************************************************\
    * Fun-02 Sec-02 Sub-03:
    *  - load data in vector for sequence entry
    \*****************************************************/

    *posUL += findFirstTrueI8(resMaskI8);
    ++(*posUL);

    if(*numReadCharUL - *posUL < vectorBytes)
    { // IF I need to have an offset
      // Check if have more file to read in
      if(*numReadCharUL < lenBuffI) return 4;

      if(outFILE != 0)
      { // IF I am printing out the entry
        *(buffCStr + *posUL - 1) = '\0';
        fprintf(outFILE, "%s\n", buffCStr + startPosUL);
        //fwrite(
        //  buffCStr + startPosUL,
        //  sizeof(char),
        //  *numReadCharUL - startPosUL,
        //  outFILE
        //);

        startPosUL = 0;
      } // IF I am printing out the entry

      // Get back to char after '\n'
      fseek(
        fqFILE,
        -1 * (*numReadCharUL - *posUL),
        SEEK_CUR
      );

      *posUL = 0;

      *numReadCharUL =
        fread(buffCStr, sizeof(char), lenBuffI, fqFILE);
    } // IF I need to have an offset

    mmLoadUI8(lineVect, (buffCStr + *posUL));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-02 Sec-03:
    ^  - Move past & find number new lines in sequence
    ^  o fun-02 sec-03 sub-01:
    ^    - Get to the end of the sequence line
    ^  o fun-02 sec-03 sub-02:
    ^    - Finalize the new line count
    ^  o fun-02 sec-03 sub-03:
    ^    - Load data into vector for the spacer (+) entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-02 Sec-03 Sub-01:
    *  - Get to the end of the sequence line
    \*****************************************************/

    mmCmpEqI8(resMaskI8, lineVect, plusVect);

    while(checkVectMaskTrueI8(resMaskI8) == 0)
    { // While I am on the sequence entry
        if(*posUL + (vectorBytes << 1) >= *numReadCharUL)
        { // If I need to read in more input
            if(*numReadCharUL < lenBuffI)
                return 4; // Fastq missing q-score entry

            // Check if going to miss any new lines
            mmCmpEqI8(resMaskI8, lineVect, newLineVect);

            if(checkVectMaskTrueI8(resMaskI8))
              numSeqLinesL += findNumTruesI8(resMaskI8);

            if(outFILE != 0)
            { // If I am writing output to a file
              *(buffCStr + *posUL + vectorBytes) = '\0';
              fprintf(outFILE,"%s",buffCStr+startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *numReadCharUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            if(*posUL + vectorBytes < *numReadCharUL)
              /*Move to charaters that were not loaded*/
              fseek(
                fqFILE,
                -1 * (*numReadCharUL-*posUL-vectorBytes),
                SEEK_CUR
              );

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            goto seqEntryLoad;
        } // If I need to read in more input

        mmCmpEqI8(resMaskI8, lineVect, newLineVect);

        // Count the number of new lines in the sequence
        // entry
        if(checkVectMaskTrueI8(resMaskI8))
           numSeqLinesL += findNumTruesI8(resMaskI8);

        *posUL += vectorBytes;

        seqEntryLoad:

        mmLoadUI8(lineVect, (buffCStr + *posUL));
        mmCmpEqI8(resMaskI8, lineVect, plusVect);
    } // While I am on the sequence entry

    /*****************************************************\
    * Fun-02 Sec-03 Sub-02:
    *  - Finalize the new line count
    \*****************************************************/

    mmCmpEqI8(resMaskI8, lineVect, plusVect);

    // Check to see if any new lines
    lineVect = zeroAfterFirstTrueI8(lineVect, resMaskI8);
     // Savgin resMaskI8, so I can avoid an extra
     // comparison later
    mmCmpEqI8(tmpMaskI8, lineVect, newLineVect);
       // Count the number of sequence new lines
       // I am using zeroAFterFirstVect to remove any
       // new lines after the plus

    if(checkVectMaskTrueI8(tmpMaskI8))
      numSeqLinesL += findNumTruesI8(tmpMaskI8);

    /*****************************************************\
    * Fun-02 Sec-03 Sub-03:
    *  - Load data into vector for the spacer (+) entry
    \*****************************************************/

    *posUL += findFirstTrueI8(resMaskI8);
    ++(*posUL);

    if(*numReadCharUL - *posUL < vectorBytes)
    { // IF I need to have an offset
      // Check if have more file to read in
      if(*numReadCharUL < lenBuffI) return 4;

      if(outFILE != 0)
      { // IF I am printing out the entry
        *(buffCStr + *posUL - 1) = '\0'; // is a '+'
        fprintf(outFILE, "%s+", buffCStr + startPosUL);

        //fwrite(
        //  buffCStr + startPosUL,
        //  sizeof(char),
        //  *numReadCharUL - startPosUL,
        //  outFILE
        //);

        startPosUL = 0;
      } // IF I am printing out the entry

      /*Move to charaters that were not loaded*/
      fseek(
        fqFILE,
        -1 * (*numReadCharUL - *posUL),
        SEEK_CUR
      );

      *posUL = 0;

      *numReadCharUL =
        fread(buffCStr, sizeof(char), lenBuffI, fqFILE);
    } // IF I need to have an offset

    mmLoadUI8(lineVect, (buffCStr + *posUL));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-02 Sec-04:
    ^  - Move past spacer entry
    ^  o fun-02 sec-04 sub-01:
    ^    - Move past the spacer entry
    ^  o fun-02 sec-04 sub-02:
    ^    - Load data into vector for the q-score entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-02 Sec-04 Sub-01:
    *  - Move past the spacer entry
    \*****************************************************/
 
    mmCmpEqI8(resMaskI8, lineVect, newLineVect);

    while(checkVectMaskTrueI8(resMaskI8) == 0)
    { /*While on the spacer entry*/

        if(*posUL + (vectorBytes << 1) >= *numReadCharUL)
        { /*If ran out of buffer*/

            if(*numReadCharUL < lenBuffI)
                return 4;         // invalid fastq entry

            if(outFILE != 0)
            { // If I am writing output to a file
              *(buffCStr + *posUL + vectorBytes) = '\0';
              fprintf(outFILE,"%s",buffCStr+startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *numReadCharUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            if(*posUL + vectorBytes < *numReadCharUL)
              /*Move to charaters that were not loaded*/
              fseek(
                fqFILE,
                -1 * (*numReadCharUL-*posUL-vectorBytes),
                SEEK_CUR
              );

            *numReadCharUL =
              fread(buffCStr,
                sizeof(char),
                lenBuffI,
                fqFILE
            ); /*Read in more of the file*/

            // Make sure buffer is a c-string
            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            goto plusEntryLoad;
        } /*If ran out of buffer*/

        *posUL += vectorBytes;

        plusEntryLoad:

        mmLoadUI8(lineVect, (buffCStr + *posUL));
        mmCmpEqI8(resMaskI8, lineVect, newLineVect);
    } /*While on the spacer entry*/

    /*****************************************************\
    * Fun-02 Sec-04 Sub-02:
    *  - Load data into vector for the q-score entry
    \*****************************************************/

    *posUL += findFirstTrueI8(resMaskI8);
    ++(*posUL);

    if(*numReadCharUL - *posUL < vectorBytes)
    { // IF I need to have an offset
      // Check if have more file to read in
      if(*numReadCharUL < lenBuffI) return 4;

      if(outFILE != 0)
      { // IF I am printing out the entry
        *(buffCStr + *posUL - 1) = '\0';
        fprintf(outFILE, "%s\n", buffCStr + startPosUL);

        //fwrite(
        //  buffCStr + startPosUL,
        //  sizeof(char),
        //  *numReadCharUL - startPosUL,
        //  outFILE
        //);

        startPosUL = 0;
      } // IF I am printing out the entry

      // Get back to char after '\n'
      fseek(
        fqFILE,
        -1 * (*numReadCharUL - *posUL),
        SEEK_CUR
      );

      *posUL = 0;

      *numReadCharUL =
        fread(buffCStr, sizeof(char), lenBuffI, fqFILE);
    } // IF I need to have an offset

    mmLoadUI8(lineVect, (buffCStr + *posUL));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-02 Sec-05:
    ^  - Move past q-score entry
    ^  - Q-score etnry should have the same number of lines
    ^    as the sequence entry
    ^  o fun-02 sec-05 sub-01:
    ^    - Find the start of the next read
    ^  o fun-02 sec-05 sub-02:
    ^    - Adjust positoin for the next read
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-02 Sec-05 Sub-01:
    *  - Find the start of the next read
    \*****************************************************/

    while(numSeqLinesL > 0)
    { /*While have q-score entry lines to print out*/

        if(*posUL + (vectorBytes << 1) >= *numReadCharUL)
        { /*If ran out of buffer*/

          // Check if we hit a new line
          mmCmpEqI8(resMaskI8, lineVect, newLineVect);

           if(checkVectMaskTrueI8(resMaskI8))
           { // If I had new lines
             numSeqLinesL -= findNumTruesI8(resMaskI8);
             if(numSeqLinesL <= 0) goto cleanUp;
           } // If I had new lines

            if(outFILE != 0)
            { // If I am writing output to a file
              *(buffCStr + *posUL + vectorBytes) = '\0';
              fprintf(outFILE,"%s",buffCStr+startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *numReadCharUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            if(*numReadCharUL < lenBuffI)
            { /*If at end of the file*/
                if(numSeqLinesL > 1) return 4; 
                else if(outFILE != 0) return 1; 
            } /*If at end of the file*/

            if(*posUL + vectorBytes < *numReadCharUL)
              /*Move to charaters that were not loaded*/
              fseek(
                fqFILE,
                -1 * (*numReadCharUL-*posUL-vectorBytes),
                SEEK_CUR
              );

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            *posUL = 0;

            goto qEntryLoad;
        } // If I need to read in more input

        mmCmpEqI8(resMaskI8, lineVect, newLineVect);

        // Count the number of new lines remaing in the
        // q-score entry

        if(checkVectMaskTrueI8(resMaskI8))
        { // If I have a new line
          numSeqLinesL -= findNumTruesI8(resMaskI8);
          if(numSeqLinesL <= 0) goto cleanUp;
        } // If I have a new line

        *posUL += vectorBytes;

        qEntryLoad:

        mmLoadUI8(lineVect, (buffCStr + *posUL));
    } /*While have q-score entry lines to print out*/

    /*****************************************************\
    * Fun-02 Sec-05 Sub-02:
    *  - Adjust positoin for the next read
    \*****************************************************/

    cleanUp:

    mmCmpEqI8(resMaskI8, lineVect, newLineVect);
    *posUL += findFirstTrueI8(resMaskI8);
    ++(*posUL);

    if(outFILE != 0)
    { // If printin out the fastq entry
      *(buffCStr + *posUL - 1) = '\0';
      fprintf(outFILE, "%s\n", buffCStr + startPosUL);

      //fwrite(
      //    buffCStr + startPosUL,
      //    sizeof(char),
      //    *posUL - startPosUL,
      //    outFILE
      //);
    } // If printin out the fastq entry

    if(*numReadCharUL - *posUL < vectorBytes)
    { // IF I need to have an offset
      // Check if have more file to read in
      if(*numReadCharUL < lenBuffI) return 1;

      // Get back to char after '\n'
      fseek(
        fqFILE,
        -1 * (*numReadCharUL - *posUL),
        SEEK_CUR
      );

      *posUL = 0;

      *numReadCharUL =
        fread(buffCStr, sizeof(char), lenBuffI, fqFILE);
    } // IF I need to have an offset

    return 0;
} /*moveToNextFastqEntry*/

// Non-vector functions
#else

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: moveToNextFastqEntry
   '  - Moves the the next fastq entry and prints to the
   '    output file if an output file was provided.
   '  - This is the vector free version of
   '    moveToNextFastqEntry
   '  o fun-04 sec-01:
   '    -  Variable declerations
   '  o fun-04 sec-02:
   '    - Move past header
   '  o fun-04 sec-03:
   '    - Find number newlines in sequence & move to spacer
   '  o fun-04 sec-04:
   '    - Move past spacer
   '  o fun-04 Sec-05:
   '    - Move past Q-score entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned long startPosUL = *posUL;
    long numSeqLinesL = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-02:
    ^  - Move past header
    ^  o fun-04 sec-02 sub-01:
    ^    - Initial set up / ensure have data in buffer
    ^  o fun-04 sec-02 sub-02:
    ^    - Move past the header
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-04 Sec-02 Sub-01:
    ^  - Initial set up / ensure have data in buffer
    \*****************************************************/

    if(*posUL >= *numReadCharUL)
    { // IF I need to get more input
      if(*numReadCharUL < lenBuffI) return 1; // EOF

      /*Read in more of the file*/
      *numReadCharUL =
        fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

      *posUL = 0;
    } // IF I need to get more input

    if(*posUL == 0) startPosUL = 0;

    /*****************************************************\
    * Fun-04 Sec-02 Sub-02:
    *  - Move past the header
    \*****************************************************/

    while(*(buffCStr + *posUL) != '\n')
    { /*While on header, move to sequence line*/

        if(*posUL >= *numReadCharUL)
        { /*If ran out of buffer */
            if(*numReadCharUL < lenBuffI)
                return 4;      // incomplete fastq entry

            if(outFILE != 0)
            { // If I am writing output to a file
              fprintf(outFILE,"%s",buffCStr + startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *posUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            // Keep as a c-string
            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            continue;
        } /*If ran out of buffer*/

        ++(*posUL);
    } /*While on header, move to sequence line*/

    ++(*posUL); // Get off the new line

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-03:
    ^  - Move past & find number new lines in sequence
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*(buffCStr + *posUL) != '+')
    { // While I am on the sequence entry
        if(*posUL >= *numReadCharUL)
        { // If I need to read in more input
            if(*numReadCharUL < lenBuffI)
                return 4; // Fastq missing q-score entry

            if(outFILE != 0)
            { // If I am writing output to a file
              fprintf(outFILE,"%s",buffCStr + startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *posUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            continue;
        } // If I need to read in more input

        // Count the number of new lines in the sequence
        // entry
        if(*(buffCStr + *posUL) == '\n') ++numSeqLinesL;

        ++(*posUL);
    } // While I am on the sequence entry

    ++(*posUL); // Get off the '+'

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-04:
    ^  - Move past spacer entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    while(*(buffCStr + *posUL) != '\n')
    { /*While on the spacer entry*/

        if(*posUL >= *numReadCharUL)
        { /*If ran out of buffer*/

            if(*numReadCharUL < lenBuffI)
                return 4;         // invalid fastq entry

            if(outFILE != 0)
            { // If I am writing output to a file
              fprintf(outFILE,"%s",buffCStr + startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *posUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            *numReadCharUL =
              fread(buffCStr,
                sizeof(char),
                lenBuffI,
                fqFILE
            ); /*Read in more of the file*/

            // Make sure buffer is a c-string
            *(buffCStr + *numReadCharUL) = '\0';
            *posUL = 0;

            continue;
        } /*If ran out of buffer*/

        ++(*posUL);
    } /*While on the spacer entry*/

    ++(*posUL); // Get off the new line

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-04 Sec-05:
    ^  - Move past q-score entry
    ^  - Q-score etnry should have the same number of lines
    ^    as the sequence entry
    ^  o fun-04 sec-05 sub-01:
    ^    - Find the start of the next read
    ^  o fun-04 sec-05 sub-02:
    ^    - Adjust positoin for the next read
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************************************************\
    * Fun-04 Sec-05 Sub-01:
    *  - Find the start of the next read
    \*****************************************************/

    while(numSeqLinesL > 0)
    { /*While have q-score entry lines to print out*/

        if(*posUL >= *numReadCharUL)
        { /*If ran out of buffer*/
            if(outFILE != 0)
            { // If I am writing output to a file
              fprintf(outFILE,"%s",buffCStr + startPosUL);
              //fwrite(
              //  buffCStr + startPosUL,
              //  sizeof(char),
              //  *posUL - startPosUL,
              //  outFILE
              //);

              startPosUL = 0;
            } // If I am writing output to a file

            if(*numReadCharUL < lenBuffI)
            { /*If at end of the file*/
                if(numSeqLinesL > 1) return 4; 
                else return 1;  // finshed the fastq file
            } /*If at end of the file*/

            /*Read in more of the file*/
            *numReadCharUL =
              fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

            *posUL = 0;

            continue;
        } // If I need to read in more input

        // Count the number of new lines remaing in the
        // q-score entry

        if(*(buffCStr + *posUL) == '\n')
          --numSeqLinesL;

        ++(*posUL);
    } /*While have q-score entry lines to print out*/

    /*****************************************************\
    * Fun-04 Sec-05 Sub-02:
    *  - Adjust positoin for the next read
    \*****************************************************/

    cleanUp:

    if(outFILE != 0)
    { // If printin out the fastq entry
      *(buffCStr + *posUL - 1) = '\0';
      fprintf(outFILE,"%s\n",buffCStr + startPosUL);

      //fwrite(
      //  buffCStr + startPosUL,
      //  sizeof(char),
      //  *posUL - startPosUL,
      //  outFILE
      //);
    } // If printin out the fastq entry

    if(*posUL >= *numReadCharUL)
    { // IF I need to get more input
      if(*numReadCharUL < lenBuffI) return 1; // EOF

      /*Read in more of the file*/
      *numReadCharUL =
        fread(buffCStr,sizeof(char),lenBuffI,fqFILE);

      *posUL = 0;
    } // IF I need to get more input

    return 0;
} /*moveToNextFastqEntry*/

#endif
