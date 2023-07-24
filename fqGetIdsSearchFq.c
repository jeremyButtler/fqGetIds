/*##############################################################################
# Name: fastqGrepSearchFastq
# Use:
#    Extracts target reads from fastq file
# Requires:
#    fastqGrepAVLTree
#    fastqGrepStructs (called by fastqGrepAVLTree)
##############################################################################*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fastqGrepSearchFastq TOC:
#    fun-1: fastqExtractTree: 
#    fun-2: buildAvlTree: Builds an balanced tree of read id's from filter file
#    fun-3: extractReads: Extract reads from fastq file with tree/hash
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include "fqGetIdsSearchFq.h"

/*---------------------------------------------------------------------\
| Output:
|   - Stdout: Prints out kept reads
|   - Returns:
|     o 0 for no errors
|     o 2 if invalid filter file
|     o 4 if invalid input fastq file
|     o 8 if could not open the output file
|     o 16 if both filter and fastq file coming from stdin
|     o 64 if had a memory error
\---------------------------------------------------------------------*/
uint8_t fastqExtract(
    char *filtPathCStr,        /*Path to file with read ids to extract*/
    char *fqPathCStr,          /*Path to fastq file to extract reads*/
    char *outPathCStr,         /*Path to fastq file to to write reads*/
    uint8_t sizeReadStackUChar,/*Number of elements to use in stack*/
    uint32_t lenBuffUI,     /*Size of buffer to read input with*/
    uint8_t hashSearchChar,    /*1: do hash search, 0: do Tree search*/
    uint8_t printReverseChar   /*1: Keep reads in filter file
                                 0: ingore reads in filter file*/
) /*Searches and extracts reads from a fastq file using read id's*/
{ /*fastqExtract*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-1 TOC: Pull out reads of interest
    '    fun-1 sec-1: Variable declerations
    '    fun-1 sec-2: Check if files exist
    '    fun-1 sec-3: Build the tree for reads
    '    fun-1 sec-4: Call tree or hash table function to search file
    '    fun-1 sec-5: Handle errors, clean up, & exit
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char buffCStr[lenBuffUI + vectorBytes + 1];  /*Holds the input line*/

    /*Majic number for kunth multiplicative hashing*/
    unsigned long majicNumULng = 0;
    uint8_t digPerKeyUChar = 0;  /*Number digitis used per key in hash*/

    /*Assume hash tabel failed, makes AVL error out when fail*/
    uint8_t hashFailedBool = 1;
    uint64_t hashSizeULng = 0;       /*Will hold Size of hash table*/

    uint64_t fastqErrULng = 0;       /*Tells if error in fastq entry*/

    /*Stack to use for searching tree. (depth = 73 = 10^18 nodes)*/
    struct readNodeStack readStack[sizeReadStackUChar + 2];
    struct readInfo *readTree = 0;
    struct readInfo **hashTbl = 0;

    FILE *filtFILE = 0;          /*to file with ID's to search for*/
    FILE *fqFILE = 0; /*fastq file to search*/
    FILE *outFILE = 0;           /*File to write extracted reads to*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-2: Check if files exist
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(filtPathCStr == 0 && fqPathCStr == 0)
        return 16;

    if(filtPathCStr == 0)
         filtFILE = stdin;
    else
        filtFILE = fopen(filtPathCStr, "r");

    if(filtFILE == 0)
        return 2;      /*No input file with read ids*/

    if(fqPathCStr == 0)
        fqFILE = stdin;
    else
        fqFILE = fopen(fqPathCStr, "r");

    if(fqFILE == 0)
    { /*If could not open the fastq file*/
        fclose(filtFILE);
        return 4;      /*No input file with reads to extract*/
    } /*If could not open the fastq file*/

    if(outPathCStr == 0) 
        outFILE = stdout;
    else
        outFILE = fopen(outPathCStr, "w");

    if(outFILE == 0)
    { /*if I could not open the output file*/
        fclose(filtFILE);
        fclose(fqFILE);
    } /*if I could not open the output file*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-1 Sec-3: Build the tree or hash table for reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Make sure start & end of my stacks are marked*/
    readStack[0].readNode = 0;
    readStack[sizeReadStackUChar + 1].readNode = 0;

    if(hashSearchChar == 0)
    { /*If just using the avl tree for searching*/
        readTree =
            buildAvlTree(
                filtFILE,    /*File with target read ids*/
                readStack,   /*Stack for searching trees*/
                buffCStr,  /*Buffer to hold one line from file*/
                lenBuffUI /*Size of buffer*/
        ); /*Build the tree of reads to search*/
    } /*If just using the avl tree for searching*/

    else
    { /*Else I am searching using a hash function*/
        hashTbl =
            makeReadHash(
                filtFILE,     /*Name of file with read ids to hash*/
                buffCStr,     /*Buffer to hold one line from file*/
                lenBuffUI,   /*Size of buffer to hold each line of the file*/
                readStack,      /*Stack for searching trees in hash table*/
                &hashSizeULng,  /*Will hold Size of hash table*/
                &digPerKeyUChar, /*Number digitis used per key in hash*/
                &majicNumULng,   /*Will hold the majic number*/
                &hashFailedBool /*Holds if manged to make hash table*/
        ); /*Build the hash table*/
    } /*Else I am searching using a hash function*/

    fclose(filtFILE); /*No longer need open*/

    if(readTree == 0 && hashFailedBool == 1) /*Calloc errored out*/
        return 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-3: Call the tree or hash table function to search the file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fastqErrULng =
        extractReads(
            fqFILE,         /*to fastq file to search through*/
            outFILE,           /*File to write extracted reads to*/
            buffCStr,        /*Buffer to hold input from fastq file*/
            lenBuffUI,      /*Size of buffCStr*/
            majicNumULng,      /*Holds majick number for kunths hash*/
            digPerKeyUChar,    /*Digits needed to get a key*/
            &printReverseChar,
            readTree,          /*AVL tree to search if hashTbl == 0*/
            hashTbl            /*hash table to search*/
    );

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-1 Sec-4: Handle errors, clean up, & exit
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    fclose(fqFILE); /*No longer need open*/
    fclose(outFILE); /*No longer need open*/

    /*Check if freeing tree or hash table with tree*/
    if(hashSearchChar == 0)
        freeReadTree(&readTree, readStack);
    else
    { /*If used hashing, free the hashing variables*/
        freeHashTbl(&hashTbl, &hashSizeULng, readStack);
    } /*If used hashing, free the hashing variables*/

    if(fastqErrULng & 4) return fastqErrULng;
    else if(fastqErrULng & 64) return 64;

    return 0; /*Success*/
} /*fastqExtract*/

/*##############################################################################
# Output:
#    Returns: balanced readInfo tree with all read id's in filtFILE
#    Returns: 0 if calloc errored out
##############################################################################*/
struct readInfo * buildAvlTree(
    FILE *filtFILE,                 /*File with read ids to keep or ignore*/
    struct readNodeStack *readStack,  /*Stack to use in building AVL tree*/
    char *buffCStr,        /*Buffer to hold one line from file*/
    uint32_t lenBuffUI         /*Size of buffer to read each line*/
) /*Builds a readInfo tree with read id's in filtFILE*/
{ /*buildAvlTree function*/

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun-2 TOC:
    '    fun-2 Sec-1: Variable declerations
    '    fun-2 Sec-2: Read id's in filter file and build tree
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-1: Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    char *tmpIdCStr = 0;
    unsigned char maxHexChar = 1;  /*Size of bigNum array*/
    uint64_t lenInUL = lenBuffUI; /*length of read id*/

    struct readInfo *readTree = 0;
    struct readInfo *lastRead = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun-2 Sec-2: Read id's in filter file and build tree
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Intalize read in loop*/
   *(buffCStr + lenBuffUI) = '\0';
   tmpIdCStr = buffCStr + lenBuffUI; /*Force initial read in*/
   lenInUL = lenBuffUI; /*Make sure read in first read*/


   do { /*While ids to read in*/
       lastRead = 
           cnvtIdToBigNum(
               buffCStr,
               lenBuffUI,
               &tmpIdCStr,
               &lenInUL,
               &maxHexChar,
               filtFILE
       ); /*Read in id and convert to big number*/

       if(lastRead == 0)
       { /*If was a falied read*/
           if(lenInUL == 0)
           { /*If was a memory allocation error (message already printed)*/
               freeReadTree(&readTree, readStack);
               return 0;
           } /*If was a memory allocation error (message already printed)*/

           break; /*end of file*/
       } /*If was a falied read*/

       if(insertNodeIntoReadTree(lastRead, &readTree, readStack) == 0)
       { /*If id is in tree, need to free*/
           freeReadInfoStruct(&lastRead); /*If id already in tree*/
           lastRead = readTree;           /*Prevent loop ending early*/
       } /*If id is in tree, need to free*/

        while(*tmpIdCStr != '\n')
        { /*While not on next entry*/
            tmpIdCStr++;

            if(*tmpIdCStr == '\0')
            { /*If at the end of the buffer, but not at start of read*/
                if(lenInUL < lenBuffUI)
                    break; /*At end of file*/

                lenInUL = fread(
                                   buffCStr,
                                   sizeof(char),
                                   lenBuffUI,
                                   filtFILE
                ); /*Read in more of the file*/

                *(buffCStr + lenInUL) = '\0';/*make sure a c-string*/
                tmpIdCStr = buffCStr;
            } /*If at the end of the buffer, but not at start of read*/
        } /*While not on next entry*/
   } while(lastRead != 0);/*ids to read in*/

   return readTree;
} /*buildAvlTree function*/

/*--------------------------------------------------------\
| Output:
|  - prints 
|    o Kept reads to outFILE
|  - Returns:
|    o 0 for no errors
|    o 4 for invalid fastq file
|    o 64 for memory allocation error
\--------------------------------------------------------*/
uint8_t extractReads(
    FILE *fqFILE,            /*fastq file to search through*/
    FILE *outFILE,              /*File to write extracted reads to*/
    char *buffCStr,        /*Buffer to hold input from fastq file*/
    uint32_t lenBuffUI,      /*Size of buffCStr*/
    unsigned long majicNumULng, /*Holds majick number for kunths hash*/
    uint8_t digPerKeyUChar,     /*Digits needed to get a key*/
    uint8_t *printNonMatchBool, /*1: print non-match, 0: print match*/
    struct readInfo *readTree,  /*For AVL search (hashTbl == 0)*/
    struct readInfo **hashTbl   /*Hash table to search for ids in*/
) /*Extract target reads from fastq file with hash table or tree*/
{ /*extractReadsInHash*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 TOC:
    #    fun-3 sec-1: Variable declerations
    #    fun-3 sec-2: Extract target reads from fastq file
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-1: Variable declerations
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uint8_t eofC = 0; // Tells if error or EOF

    // I need a blank number to initalize the big number
    // structuer
    char dummyConvertUChar = '0';

    // Maximum number of limbs to have in big number struct
    unsigned int maxLimbsUI = 1024;
      // This is 16 limbs for 64 bit integers

    unsigned long lenInUL = 0; // Number of chars in buffer

    unsigned long posUL = 0;   // Position at in buffer
    unsigned long startUL = 0; // Start of fastq entry

    // For hash table
    struct bigNum *idBigNum = 0;

    // For searching read id tree (part of hash table)
    struct readInfo *lastRead = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Fun-3 Sec-2: Extract target reads from fastq file
    #    fun-3 sec-3 sub-1: Get the read name of a single fastq entry
    #    fun-3 sec-3 sub-2: Determine if read is in tree
    #    fun-3 sec-3 sub-3: Decide if should keep read
    #    fun-3 sec-3 sub-4: Keeping read, print out read
    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    // Read in the first bit of input
    lenInUL =
      fread(buffCStr, sizeof(char), lenBuffUI, fqFILE);

    buffCStr[lenBuffUI] = '\0';          /*Make c-string for fread*/
    posUL = 0;

    // Set up the big number sruct
    idBigNum =
      makeBigNumStruct(&dummyConvertUChar, &maxLimbsUI);

    if(idBigNum == 0) return 64;

    while(eofC == 0)
    { /*While there are lines in the file*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-1: Get the read name of a single fastq entry
        ***********************************************************************/

        eofC =
            parseFastqHeader(
                buffCStr, // Buffer to hold file input
                &posUL,   // Positon working on in buffer
                //&startUL, // Reset if read in new buffer
                &lenInUL,
                lenBuffUI,
                idBigNum,
                fqFILE
        ); /*Get read name from file*/

        if(startUL > posUL) startUL = 0;

         if(eofC)
         { // If at end of file
             if(eofC & 1) eofC = 0;

             freeBigNumStruct(&idBigNum);
             return eofC; // Will be 4 for invalid fastq
                          // and 1 for complete fastq
         } // If at end of file

        /***********************************************************************
        # Fun-3 Sec-3 Sub-2: Determine if read is in tree
        ***********************************************************************/

        // Check if using hash table/AVL tree or just AVL
        if(hashTbl == 0)
            lastRead = searchTree(idBigNum, readTree);
        else 
        { /*Else doing a hash table / AVL tree*/
            lastRead = 
                findReadInHashTbl(
                    idBigNum,
                    &majicNumULng,
                    &digPerKeyUChar,
                    hashTbl
            );  /*See if read id is in the hash table*/
        } /*Else doing a hash table / AVL tree*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-3: Decide if should keep read
        ***********************************************************************/

        /*Check if should print read (!! converts address to 1)
          read is match, priting matches only = 1 ^ 0 = 1
          read is match, priting non-matches only = 1 ^ 1 = 0
          read is not match, priting matches only = 0 ^ 0 = 0
          read is not match, priting non-matches only = 0 ^ 1 = 1
        */
        if(((!!lastRead) ^ *printNonMatchBool) == 0)
        { /*If is a read I am not printing out*/
             eofC =
                 moveToNextFastqEntry(
                     buffCStr,
                     &posUL,   // Points to End of read id
                     lenBuffUI,
                     &lenInUL,   /*Holds how many characters fread got*/
                     fqFILE,
                     0         // Do not print anything
             ); /*Move to next fastq entry*/


             if(eofC)
             { // If at end of file
                 if(eofC & 1) eofC = 0;

                 freeBigNumStruct(&idBigNum);
                 return eofC; // 4 for invalid fastq
                              // and 1 for complete fastq
             } // If at end of file

            startUL = posUL;
            continue;
         } /*If is a read I am not printing out*/

        /***********************************************************************
        # Fun-3 Sec-3 Sub-4: Keeping read, print out read
        ***********************************************************************/

         eofC =
             moveToNextFastqEntry(
                 buffCStr,
                 &startUL,       // Start of read id
                 lenBuffUI,
                 &lenInUL,
                 fqFILE,
                 outFILE
         ); /*Print read & move to next read*/

         posUL = startUL;

         if(eofC)
         { // If at end of file
           if(eofC & 1) eofC = 0;

           freeBigNumStruct(&idBigNum);
           return eofC; // 4 for invalid fastq
                        // and 1 for complete fastq
         } // If at end of file
    } /*While there are lines in the file*/

    freeBigNumStruct(&idBigNum);
    return 0; /*No errors*/
} /*extractReadsInHash*/
