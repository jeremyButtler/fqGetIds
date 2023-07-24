/*#########################################################
# Name: vectorWrap
# Use:
#  - Is a wraper for SSE2 equivlent instruction sets
#    for SSE2, AVX2, AVX512, and NEON SMID instructions.
#  - All vector functions wrapped by macros to ensure all
# - This is not designed to replace better options, such
#   as SIMDe, but it is designed to provide single
#   functions that can be complied for all supported
#   instruction sets.
#    code is inlined
#  - WARNING: Not all comparisons returen the same result.
#    o This can be avoided for bit counts by calling
#      fixIXCmpMaskCnt, which corrects the bit count.
# Compiling for gcc:
#  - sse: -DSSE
#  - avx2: -DAVX2 -mavx2
#  - avx512: -DAVX512 -mavx512BW
#    - Do not use -mavx512f.
#  - neon: -mfpu=neon-vfpv4 or -mfpu=neon
# Supporting functions:
#  - Are here to to provide support for the vectors or do
#    things that would be difficult for vecotors. You can 
#    find these in sec-04.
#    o These are mostly for int8 vectors.
#  - These functions may be slower than doing scalar
#    operations instead.
# Datatype:
#  - This is the worst world between ARM and AVX512.
#    Inorder to support all datatypes I need to support the
#    datatypes unique to NEON and AVX512.
#  - vectIX is an alias for the vector data type. The X in
#    vectIX is the size of the datatype working on with the
#    vecotor. For example a vectors of int8_t's would be
#    vectI8.
#    o Currnetly datatypes for vectI8, vectI16, and
#      vectI32.
#    o 64 bit numbers are not supported, due to intel
#      having no comaparisons for SSE2.
#  - For comparisions their are three additional data types
#    o mmaskI8
#       - for __mmask64 data type (8 bit comparisons)
#      - This is returned for all 8 bit comparisons
#    o mmaskI16
#      - for __mmask32 data type (16 bit comparisons)
#      - This is returned for all 16 bit comparisons
#      - 16 bit comparisoins will likely behave differntly
#        in AVX512 than AVX2 or SSE.
#        - This is due to a change in intels return type
#    o mmaskI32
#      - for __mmask16 data type (32 bit comparisons)
#      - This is returned for all 32 bit comparisons
#      - 32 bit comparisoins will likely behave differntly
#        in AVX512 than AVX2 or SSE.
#        - This is due to a change in intels return type
#  - One warning about macros:
#    Sending in a pointer by "pointer + newPosition" does
#    not work properly if parentheses are not used. In
#    this case do (pionter + newPosition)
# Definitions:
#   - vectorBits: Is the number of bits in the vector
#   - vectorBytes: Is the number of bytes a vector holds
# Wrapper functions:
#  - For each function there are 8 bit, 16 bit, and 23 bit
#    variants. This is represented by X in the function
#    namnes listed here.
#  - Functions either take ro return a mmaskIX or vactIX.
#  - Input/Output
#    o mmLoadIX
#      - Loads a an array into the vector
#      - Wrapper for _mmx_load_siX((__mxi *) array)
#    o mmLoadUIX
#      - mmLoad for unaligned data
#    o mmSetZeroIX
#      - Sets all values in a vectors with zeros
#      - Wrapper for _mmx_setzero_siX()
#      - Note: this is equivlent to mmSet1IX for NEON.
#        So, do not depend on it always being a latency
#        of 1, like for SSE/AVX.
#    o mmSet1IX
#      - Sets all values in a vectors to input number
#      - Wrapper for _mmx_set1_epiX()
#    o mmStoreIX
#      - Stores a vector into an array
#      - Use vectorBytes to determine the array size in
#        bytes
#      - Wrapper for mmx_store_siX((__mxi *) array, vector)
#      - This function casts input array, so it is only
#        datatype specific for the vector type, not the
#        array.
#  - Compaisons (returns mmaskx)
#    o fixIXCmpMaskCnt:
#      - These functions make sure that bit counts are
#        corrected for the differences in return
#        comparisons between comparisons.
#    o mmStoreIXCmpMask:
#      - This stores a mmaskIX into an unsigned long
#      - Note, because of differences between SSE, AVX512,
#        and NEON, this function will have slighty
#        different output.
#        o NEON64 returns the acutal vector full of ones
#          or zeros.
#        o NEON returns a long with multiple bits per
#          datatype (number bits in datatype / 2).
#        o SSE and AVX2 returns an long (techincally int)
#          with 1 bit for every 8 bits
#          (number of bits in datatype / 8).
#        o AVX512 returns a mask with 1 bit per datatype
#       - You can correct counts and a bit index by calling
#         fixIXCmpMaxkCnt
#    o mmCmpGtIX
#      - Do a greater then comparison on vectors
#      - Returns an mmaskIX
#    o mmCmpLtIX
#      - Do a lesser then comparison on vectors
#      - Returns an mmaskIX
#    o mmCmpEqIX
#      - Do an is equal comparison on vectors
#      - Returns an mmaskIX
#  - Logical
#    o mmAndNotIX
#      - Not the first vector and then do an & (and) on the
#        first and second vector
#      - Wrapper for mmX_andnot_siX(vect1, vect2)
#    o mmAndIX
#      - Perform an & (and) opteration on both vectors
#      - Wrapper for mmX_and_siX(vect1, vect2)
#    o mmOrIX
#      - Perform an | (or) opteration on both vectors
#      - Wrapper for mmX_or_siX(vect1, vect2)
#    o mmXOrIX
#      - Perform an ^ (xor) opteration on both vectors
#      - Wrapper for mmX_or_siX(vect1, vect2)
#  - Bit mainpulation
#    o mmShiftRightIX
#      - shift a vector right by input bytes
#    o mmShiftLeftIX
#      - shift vector left by input bytes
#  - Math
#    o mmAddIX
#      - adds two vectors together
#    o mmAddSatIX
#      - Saturation add two vectors
#    o mmSubIX
#      - subtact two vectors
#    o mmSubSatIX
#      - Satuturation Subtact two vectors
#  - Vector conversons (casting)
#    o cnvtIXToIY
#      - Converts one vectIX type to a vectIY type.
#      - EX: convert vectI8 to vectI16
# Libraries:
# C Standard Libraries:
#   - <limits.h>
#   - <stdint.h>
#   - <immintrin.h>
#   - <arm_neon.h>
# TODO:
#   - Test arm_neon instructions
#########################################################*/

#include "vectorWrap.h"

// Disable functions if not using vectors
// This is to avoide makfile errors when not using vectors
#if defined AVX512 || AVX2 || SSE || NEON || NEON64

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'  o fun-01 makeVectBufferMalloc:
'     - Make an aligned buffer for a vector
'  o fun-02 makeVectBufferMalloc:
'     - Make an aligned buffer with zeros for a vector
'  o fun-03 makeVectBufferMalloc:j
'     - Resize a buffer making sure it is aligned to an
'       vector using realloc
'  o fun-04 bytesToAllocate:
'    - Find the number of bytes to allocate
'  o fun-05 checkVectMaskFalse:
'    - Check to see if a mask from a vector comparison at
'      least one false value
'  o fun-06 checkVectMaskTrue:
'    - Check to see if a vector mask has at least one true
'      value (1).
'  o fun-07 findFirstFalse:
'     - Finds the first zero char value in a vector mask
'  o fun-08 findLastTrue:
'    - Finds the last true value (1) in a comparison mask
'      from comparing two vectors
'  o fun-09 findFirstTrue:
'     - Finds the first true value (1) in a comparison mask
'       from comparing two vectors
'  o fun-10 findNumTrues:
'    - Find the number of set bits (1's or true's) in a
'      vector.
'  o fun-11 zeroAfterFirstTrue:
'    - Zeros all elements in the vector after the first
'      true (1) value
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Pointer to buffer the size of the vector
|    o 0 if malloc had a memeory error
\--------------------------------------------------------*/
void * makeVectBufferMalloc(
  unsigned long numBytesUL
   // Number of elements to assign to buffer
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: makeVectBufferMalloc
   '  - Make an aligned buffer for a vector
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  numBytesUL += vectorBytes - (numBytesUL % vectorBytes);
  return malloc(sizeof(char) * numBytesUL);
} // makeVectBufferMalloc

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Pointer to buffer the size of the vector
|    o 0 if calloc had a memeory error
\--------------------------------------------------------*/
void * makeVectBufferCalloc(
  unsigned long numBytesUL
   // Number of elements to assign to buffer
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: makeVectBufferMalloc
   '  - Make an aligned buffer with zeros for a vector
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  numBytesUL += vectorBytes - (numBytesUL % vectorBytes);
  return calloc(numBytesUL, sizeof(char));
} // makeVectBufferMalloc

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o Pointer to the resized buffer
|    o 0 if realloc had a memeory error
\--------------------------------------------------------*/
void * makeVectBufferRealloc(
  void *buffToResize, // Buffer to resize to new size
  unsigned long numBytesUL
   // Number of elements to assign to buffer
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: makeVectBufferMalloc
   '  - Resize a buffer making sure it is aligned to an
   '    vector using realloc
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  numBytesUL += vectorBytes - (numBytesUL % vectorBytes);
  return realloc(buffToResize, numBytesUL * sizeof(char));
} // makeVectBufferMalloc

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o unsinged long with number of bytes to use in buffer
\--------------------------------------------------------*/
unsigned long bytesToAllocate(
  unsigned long numBytesUL
   // Planned array size
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: bytesToAllocate
   '  - Find the number of bytes to allocate
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  return
    numBytesUL + vectorBytes - (numBytesUL % vectorBytes);
} // bytesToAllocate

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: If there were no zeros in the vector mask
|    o 1: If there was a zero (false) in the vector mask
\--------------------------------------------------------*/
char checkVectMaskFalseI8(
  mmaskI8 vectMask // Mask to check for zeros
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: checkVectMaskFalse
   '  - Check to see if a mask from a vector comparison at
   '    least one false value
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  // Latency of 3 for 256 and 128; 6 to 8 for 512 bit
  storeI8CmpMask(maskUX, vectMask);

  #ifdef AVX512
    if(maskUX < defMaxUI64) return 1;
  #elif AVX2
    if(maskUX < defMaxUI32) return 1;
  #elif SSE
    //maskUX = maskUX >> 16; // Remove 16 zeros
      // Turns out that the code fills up the lowest
      // postions first, so I did not need this
    if(maskUX < defMaxUI16) return 1;
  #endif

  return 0; // No zeros in the vector
} // checkVectMaskFalse

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: If there is only false values
|    o 1: If there is a true value
\--------------------------------------------------------*/
char checkVectMaskTrueI8(
  mmaskI8 vectMask // Mask to check for zeros
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: checkVectMaskTrue
   '  - Check to see if a vector mask has at least one true
   '    value (1).
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  // Latency of 3 for 256 and 128; 6 to 8 for 512 bit
  storeI8CmpMask(maskUX, vectMask);

  if(maskUX != 0) return 1;
  return 0; // No zeros in the vector
} // checkVectMaskTrue

// Need to add in logic to handle OS were 64 bits varibles
// are not present
/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o char pointer to first zero in the mask or to the
|      last base in the c-string (posCStr)
| Note:
|  - Cycles:
|    - This function calls _mmx_store, which has a latency
|      of 5 on skylake. It then checks all longs to see
|      if the first part has a zero. After that it does
|      a binary search for the first zero in the long.
| WARNING
|   - I have not tested this code yet
\--------------------------------------------------------*/
unsigned char findFirstFalseI8(
  mmaskI8 vectorMask      // vector mask to find zero in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: findFirstFalse
   '  - Finds the first zero char value in a vector mask
   '    and returns a pointer to it or to the end of the
   '    string
   '  o fun-07 sec-01:
   '    - Varaible declerations
   '  o fun-07 sec-02:
   '    - Find the position of the first zero in the long
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-07 Sec-01:
  ^  - Varaible declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned char posUC = 0;

  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-07 Sec-02:
  ^  - Find the position of the first zero in the long
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  // Not sure if I have this right, but doing backwards
  // because intel loads its vectors backwards. So the
  // start of the string is the smallest bit
  // Latency of 3 for 256 and 128; 6 to 8 for 512 bit
  storeI8CmpMask(maskUX, vectorMask);

  // If dealing with a 64 bit mask
  #if defined AVX512 || NEON || NEON64
    // ADD IFNDEF 64 BIT OS
    if(maskUX >> 32 != defMaxUI32)
    { // If all last digits set (vectors are backwards)
      maskUX = maskUX >> 32;
      posUC += 32;
    } // If all last digits set (vectors are backwards)

    else maskUX &= 0x00000000FFFFFFFF;
  #endif

  // Not needed for 128 bit vectors
  #if defined AVX2 || AVX512 ||  NEON || NEON64
    if(maskUX >> 16 != defMaxUI16)
    { // If all last digits set (vectors are backwards)
      maskUX = maskUX >> 16;
      posUC += 16;
    } // If all last digits set (vectors are backwards)

    else maskUX &= 0x0000FFFF;

  //#elif SSE
    // Fills form lowest position up, so do not need
    //maskUX = maskUX >> 16; // Lower 16 is just 0's
  #endif

  if(maskUX >> 8 != defMaxUI8)
  { // Else zero in last 1/4 of vector (first 1/4 of buff)
    maskUX = maskUX >> 8;
    posUC += 8;
  } // Else zero in last 1/4 of vector (first 1/4 of buff)

  else maskUX &= 0x00FF;

  if(maskUX >> 4 != 0xF)
  { // If the next 4 digits are all ones
    maskUX = maskUX >> 4;
    posUC += 4;
  } // If the next 4 digits are all ones

  else maskUX &= 0x0F;

  // Check To see which two bits it is
  if(maskUX >> 2 != 3)
  { // If the next 2 digits are all ones
    maskUX = maskUX >> 2;
    posUC += 2;
  } // If the next 2 digits are all ones

  else maskUX &= 3; // Keep only the 2nd ahd thrid bit

  // Check the final value
  if(maskUX >> 1 == 0) ++posUC;

  // Adust for off counts
  fixI8CmpMaskCnt(posUC, posUC);

  return posUC;
} // findFristFalse

// Need to add in logic to handle OS were 64 bits varibles
// are not present
/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o char pointer to the last true (1) in the mask or
|      to the last base in the c-string (posCStr)
| Note:
|  - Cycles:
|    - This function calls _mmx_store, which has a latency
|      of 5 on skylake. It then checks all longs to see
|      if the first part has a zero. After that it does
|      a binary search for the first zero in the long.
\--------------------------------------------------------*/
unsigned char findLastTrueI8(
  mmaskI8 vectorMask     // vector mask to find zero in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-08 TOC: findLasttTrue
   '  - Finds the last true value (1) in a comparison mask
   '    from comparing two vectors
   '  o fun-08 sec-01:
   '    - Varaible declerations
   '  o fun-08 sec-02:
   '    - Find the position of the first zero in the long
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-08 Sec-01:
  ^  - Varaible declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned char posUC = 0;

  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-08 Sec-02:
  ^  - Find the position of the first zero in the long
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  // Latency of 3 for 256 and 128; 6 to 8 for 512 bit
  storeI8CmpMask(maskUX, vectorMask);

  // Doing checks backwards because the vectors are loaded
  // backwards. The start of the string is the smallest
  // bit
  // If dealing with a 64 bit mask
  #if defined AVX512 ||  NEON || NEON64
    // ADD IFNDEF 64 BIT OS
    if(maskUX >> 32 != 0)
    { // If the position is in the higher bit
      maskUX = maskUX >> 32;
      posUC += 32;
    } // If the position is in the higher bit
  #endif

  // Not needed for 128 bit vectors
  #if defined AVX2 || AVX512 ||  NEON || NEON64
    if(maskUX >> 16 != 0)
    { // If the position is in the higher bit
      maskUX = maskUX >> 16;
      posUC += 16;
    } // If the position is in the higher bit

  //#elif SSE
    // Fills form lowest position up, so do not need
    //maskUX = maskUX >> 16; // Lower 16 is just 0's
  #endif

  if(maskUX >> 8 != 0)
  { // If I have a higher bit value
     maskUX = maskUX >> 8;
     posUC += 8;
  } // If I have a higher bit value

  if(maskUX >> 4 != 0)
  { // If the value is in the higher bit
    maskUX = maskUX >> 4;
    posUC += 4;
  } // If the value is in the higher bit

  // Check To see which two bits it is
  if(maskUX >> 2 != 0)
  { // If is in the higher bit
    maskUX = maskUX >> 2;
    posUC += 2;
  } // If is in the higher bit

  if(maskUX >> 1 != 0) ++posUC;

  // Adust for off counts
  fixI8CmpMaskCnt(posUC, posUC);

  return posUC;
} // findLastTrue

// Need to add in logic to handle OS were 64 bits varibles
// are not present
/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o char pointer to the first true (1) in the mask or
|      to the last base in the c-string (posCStr)
| Note:
|  - Cycles:
|    - This function calls _mmx_store, which has a latency
|      of 5 on skylake. It then checks all longs to see
|      if the first part has a zero. After that it does
|      a binary search for the first zero in the long.
\--------------------------------------------------------*/
unsigned char findFirstTrueI8(
  mmaskI8 vectorMask     // vector mask to find zero in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-09 TOC: findFirstTrue
   '  - Finds the first true value (1) in a comparison mask
   '    from comparing two vectors
   '  o fun-09 sec-01:
   '    - Varaible declerations
   '  o fun-09 sec-02:
   '    - Find the position of the first zero in the long
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-09 Sec-01:
  ^  - Varaible declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned char posUC = 0;

  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-09 Sec-02:
  ^  - Find the position of the first zero in the long
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  // Latency of 3 for 256 and 128; 6 to 8 for 512 bit
  storeI8CmpMask(maskUX, vectorMask);

  // Doing checks backwards because the vectors are loaded
  // backwards. The start of the string is the smallest
  // bit
  // If dealing with a 64 bit mask
  #if defined AVX512 ||  NEON || NEON64
    // ADD IFNDEF 64 BIT OS
    if(maskUX & 0x00000000FFFFFFFF)
      maskUX &= 0x00000000FFFFFFFF;

    else
    { // Else the first bit is in the uppler position
      posUC += 32;
      maskUX = maskUX >> 32;
    } // Else the first bit is in the uppler position

  #endif

  // Not needed for 128 bit vectors
  #if defined AVX2 || AVX512 ||  NEON || NEON64
    if(maskUX & 0x0000FFFF) maskUX &= 0x0000FFFF;

    else
    { // Else the first bit is in the uppler position
      posUC += 16;
      maskUX = maskUX >> 16;
    } // Else the first bit is in the uppler position

  //#elif SSE
    // Fills form lowest position up, so do not need
    //maskUX = maskUX << 16; // Lower 16 is just 0's
  #endif

  if(maskUX & 0x00FF) maskUX &= 0x00FF;

  else
  { // Else the first bit is in the uppler position
    posUC += 8;
    maskUX = maskUX >> 8;
  } // Else the first bit is in the uppler position

  if(maskUX & 0x0F) maskUX &= 0x0F;

  else
  { // Else the first bit is in the uppler position
    posUC += 4;
    maskUX = maskUX >> 4;
  } // Else the first bit is in the uppler position

  if(maskUX & 3) maskUX &= 3;

  else
  { // Else the first bit is in the uppler position
    posUC += 2;
    maskUX = maskUX >> 2;
  } // Else the first bit is in the uppler position

  if(!(maskUX & 1)) ++posUC;

  // Adust for off counts
  fixI8CmpMaskCnt(posUC, posUC);

  return posUC;
} // findFristTrue


/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The number of set bits (1's/true's) in the vector
|      mask
\--------------------------------------------------------*/
unsigned char findNumTruesI8(
  mmaskI8 vectorMask     // vector mask to find zero in
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-10 TOC: findNumTrues
   '  - Find the number of set bits (1's or true's) in a
   '    vector.
   '  - This coes is From
   '    https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
   '  - I would like to say I understand this code, however,
   '    I do not. I just know that it works
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


  #ifdef AVX512
    uint64_t maskUX;
  #elif AVX2
    uint32_t maskUX;
  #elif SSE
    uint32_t maskUX;
  #elif NEON
    uint64_t maskUX;
  #elif NEON64
    uint64_t maskUX;
  #endif

  storeI8CmpMask(maskUX, vectorMask);

  #if defined AVX2 || SSE
    maskUX = maskUX - ((maskUX >> 1) & 0x55555555);

    maskUX =
      (maskUX & 0x33333333) + ((maskUX >> 2) & 0x33333333);

    maskUX =
       (((maskUX + (maskUX >> 4)) & 0xF0F0F0F) * 0x1010101)
      >> 24; // << ((sizeof(unsigned int) - 1) * 8)
  #endif

  #if defined AVX512 || NEON || NEON64
    maskUX = maskUX - ((maskUX >> 1) & 0x5555555555555555);

    maskUX =
        (maskUX & 0x3333333333333333)
      + ((maskUX >> 2) & 0x3333333333333333);

    maskUX = 
      (
          ((maskUX + (maskUX >> 4)) & 0xF0F0F0F0F0F0F0F)
        * 0x101010101010101
      ) >> (sizeof(unsigned long) - 1) * 8;
  #endif

  fixI8CmpMaskCnt(maskUX, maskUX);
  return maskUX;
} // findNumTrues

/*--------------------------------------------------------\
| Output;
|  - Returns:
|    o vectI8 datatype will zeros for all bytes after the
|      first true byte
| Note:
|  - This is not very fast, but is the best I can do since
|    the shift left function requires a constant.
\--------------------------------------------------------*/
vectI8 zeroAfterFirstTrueI8(
  vectI8 vectToZero, // Vector to make zeros after 1st true
  mmaskI8 maskVect  // Has the one to zero after
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-11 TOC: zeroAfterFirstTrue
   '  - Zeros all elements in the vector after the first
   '    true (1) value
   ' Fun-11 Sec-02:
   '  - Mask building; checks that are only for AVX512
   ' Fun-11 Sec-03:
   '  - Mask building; checks unique to AVX2 and AVX512
   ' Fun-11 Sec-04:
   '  - Mask building; checks for SSE/AVX2/AVX512
   ' Fun-11 Sec-05:
   '  - Loading and applying the mask
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  vectI8 retVect;
  char shiftUC =
    vectorBytes - (findFirstTrueI8(maskVect)) - 1;
    // Find how many btyes till I read the first true
    // this will be the number of bytes to keep
    // -1 is to account for vectorBytes being index 1

  #ifdef AVX512
    uint64_t maskAryUL[8];
  #elif AVX2
    uint64_t maskAryUL[4];
  #elif SSE
    uint64_t maskAryUL[2];
  #elif NEON
    uint64_t maskAryUL[2];
  #elif NEON64
    uint64_t maskAryUL[1];
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-11 Sec-02:
  ^  - Mask building; checks that are only for AVX512
  ^  - Mask building is done backwards, becasue the longs
  ^    will be loaded in backwards
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  #ifdef AVX512
    if(shiftUC > 7)
      maskAryUL[7] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[7] = ULONG_MAX >> (shiftUC * 8);

      maskAryUL[6] = ULONG_MAX;
      maskAryUL[5] = ULONG_MAX;
      maskAryUL[4] = ULONG_MAX;
      maskAryUL[3] = ULONG_MAX;
      maskAryUL[2] = ULONG_MAX;
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;

      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;

    if(shiftUC > 7)
      maskAryUL[6] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[6] = ULONG_MAX >> (shiftUC * 8);

      maskAryUL[5] = ULONG_MAX;
      maskAryUL[4] = ULONG_MAX;
      maskAryUL[3] = ULONG_MAX;
      maskAryUL[2] = ULONG_MAX;
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;
      
      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;

    if(shiftUC > 7)
      maskAryUL[5] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[5] = ULONG_MAX >> (shiftUC * 8);
      maskAryUL[4] = ULONG_MAX;
      maskAryUL[3] = ULONG_MAX;
      maskAryUL[2] = ULONG_MAX;
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;
      
      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;

    if(shiftUC > 7)
      maskAryUL[4] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[3] = ULONG_MAX >> (shiftUC * 8);
      maskAryUL[2] = ULONG_MAX;
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;

      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-11 Sec-03:
  ^  - Mask building; checks unique to AVX2 and AVX512
  ^  - Mask building is done backwards, becasue the longs
  ^    will be loaded in backwards
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  #if defined AVX2 || AVX512
    if(shiftUC > 7)
      maskAryUL[3] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[3] = ULONG_MAX >> (shiftUC * 8);
      maskAryUL[2] = ULONG_MAX;
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;

      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;

    if(shiftUC > 7)
      maskAryUL[2] = 0;
    else
    { // Else I need to shift the enxt long
      maskAryUL[2] = ULONG_MAX >> (shiftUC * 8);
      maskAryUL[1] = ULONG_MAX;
      maskAryUL[0] = ULONG_MAX;
      
      goto loadMask;
    } // Else I need to shift the enxt long

    shiftUC -= 8;
  #endif

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-11 Sec-04:
  ^  - Mask building; checks for SSE/AVX2/AVX512
  ^  - Mask building is done backwards, becasue the longs
  ^    will be loaded in backwards
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  #if defined SSE || NEON || AVX2 || AVX512
    if(shiftUC > 7)
      maskAryUL[1] = 0;

    else
    { // Else only need to shift elemetns in the first limb
      maskAryUL[1] = ULONG_MAX >> (shiftUC * 8);
      maskAryUL[0] = ULONG_MAX;

      goto loadMask;
    } // Else only need to shift elemetns in the first limb

    shiftUC -= 8;
  #endif

  if(shiftUC > 7)
    maskAryUL[0] = 0;

  else
    maskAryUL[0] = ULONG_MAX >> (shiftUC * 8);

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-11 Sec-05:
  ^  - Loading and applying the mask
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  loadMask:

  mmLoadUI8(retVect, maskAryUL);
    // I will get some rare setfaults on AVX2 when I use
    // mmload. Some odd alignment issues here
  mmAndI8(retVect, vectToZero, retVect);

  return retVect;
} // zeroAfterFristTrue

#endif
