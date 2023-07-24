# Use

FqGetIds extracts reads from a fastq file by input read
  ids. It can be complied with vector support for
  processing the fastq file cpu's that support the NEON,
  SSE, AVX2, and maybe AVX512 instruction sets. FqGetIds
  without vector support is slower than seqkit. FqGetIds
  also requires more ram than seqkit.

I still need to benchmark to see how vector support
  compares seqkit in speed.

The reason why I am posting fqGetIds is to see if anyone
  can find a use for any elements of the code.

This program is dual licensed under the MIT (primary) or
  CC0 (alt-LICENSE). Pick the license you prefer to work
  with.

# Building fqGetIds

```
# 128 bit vector support
# should work on most if not all computers
sudo make sse
make install

# No vector (scalar) support
make
sudo make install

# 256 bit vector support
make avx2
sudo make install

# 512 bit vector support
make avx512
sudo make install

# If you get an Illegal instruction error when you run
# fqGetIds, it likely means that your cpu does not support
# the SIMD instructions you compiled fqGetIds for.
# This is most likely to happend with AVX512.

# NEON 128 bit vector support (arm cpus)
make neon
sudo make install

# NEON 64 bit vector support (128 bit is better)
make neon64
sudo make install
```

# Running fqGetIds

fqGetIds takes in a fastq file and a file of read ids as
  input. Each entry in the read id file must start with
  an readID or and @readID. The read id must be followed by
  white space, such as a space, tab, or new line. The entry
  then must end with a new line. So long as this format is
  followed your input should be ok. The code block beneath
  has some examples of what your file could look like.

```
@de83dd55-d1d3-4bd8-9de6-2f8b3930f8ec
c10ae957-ec6b-49f0-8695-7a5347afb66b
@b8b59199-b3f8-47da-82c9-93bdf6ca251a runid=...
d952424e-f1a0-4449-b56c-d8ceff123521	flag	reference	mapq	position	cigar	 ...
```

The fastq file should be a standard fastq. It can not have
  blank lines. However, it can have multiple lines for the
  sequence entry and Q-score entry, so long as the sequence
  and Q-score entry have the same number of lines.

Here are some example commands for running fqGetIds. These
  should cover all the options for fqGetIds.

```
fqGetIds -f read-ids.txt -fastq file.fastq > out.fastq

# To ingore read ids in read-ids.txt
fqGetIds -v -f read-ids.txt -fastq file.fastq -out out.fastq

# Read fastq input from stdin
cat file.fastq | fqGetIds -stdin-fastq -f read-ids.txt > out.fastq

# Read read ids to extract from stdin
cat read-ids.txt | fqGetIds -stdin-filt -fastq file.fastq -out out.fastq

# help message
fqGetIds -h | less
```

# Graphs

I will add a comparison to seqkit version 2.3 here soon.

# How fqGetIds works

fqGetIds uses a combination of a hash table and an AVL tree
  to detect if reads should be extracted or ignored. The
  first step reduces an read id to an number that can be
  hashed using Knuth's multiplicative hash. This number is
  then hashed to see if it matches any input reads. An AVL
  tree is then used to handle cases were multiple reads
  have the same hash.

In the conversion step only characters 0 to 9 and a to u
  (case is ignored) are converted to a 5 bit number (1 to
  31). Other characters of interest (v to z [case is
  ignored] and :), are converted to 0, which means that
  only their position is kept. However, this position will
  be lost if it is at the end of the last limb.

The converted read id is stored in multiple longs (each
  separate long is a limb). A sum of all limbs is found by
  casting each long limb two integer limbs (smaller than
  longs) and then adding up the integers. This sum is used
  to find the hash value and for a first comparison of two
  different read ids.

An AVL tree is used to handle collisions in the hash table.
  The first check when comparing read ids in the tree to
  a query is to make sure both reads have the same total.
  Next is the check to make sure both reads have the same
  number of limbs.Finally each limb in the read id is
  compared to see if any limbs are different. In this case
  the limbs are compared backwards (last limbs are compared
  first).

I suspect the difference in times between fqGetIds and
  seqkit are due to my AVL tree hash combination being
  slower than seqkit's hashing. I also know that Illumina
  read id's barely fill three third limbs, which results in
  the difference in memory usage between fqGetIds and
  seqkit being larger for Illumina read ids than Nanopore
  read ids.

# Vector support

FqGetIds uses a ~100kb buffer to read in sections of the
  fastq file. This reduces the number of times fqGetIds
  needs to access the original file. However, this comes
  at the cost of fqGetIds needing to scan through the
  buffer to identify the start and end of fastq files.
  To speed this step up fqGetIds can be compiled with 
  SIMD support.

Currently the NEON, SSE, and AVX2 instruction sets are
  known to work with fqGetIds. However, vectorWrap.c/h
  does support the AVX512 instructions, so it is possible
  that the AVX512 will work (I am unable to test this).

# TODO:

I am not sure how much farther I will continue this
  project and it will be a bit before I move back to
  fqGetIds.

It might be nice to try changing my AVL hash to a double
  hash. I am hoping this will get closer to seqkits times.

Also I need to add a usage guide for fqGetIds. You can find
  an old usage guide at
  [https://github.com/jeremybuttler/find--Co-infections](
  https://github.com/jeremybuttler/find--Co-infections).
  The fqGetIds entry in it is pretty brief and I am hoping
  to add a better one here.

Clean up the documentation. This program was developed over
  a couple of years, over which my documentation style has
  changed. I need to take a couple days to get the code to
  my current format (fqGetIdsFqFun is in the correct
  format).

Make fqGetIds easier to work with by using structures. This
  will only be done if it does not compromise speed, which
  it should not.

# Thanks:

My Dad for advice and support.

Danila Kutenin who wrote the blog about how to get a
  movemask function for NEON. This helped me get the NEON
  vector support set up.
  [https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon](
  https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon)

The bit hacking guide for providing a function to count
  the number of set bits.
[https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive](
 https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive)

Various internet sites, especially stack overflow, the
  intel intrinsic guide, and arms intrinsic search engine,
  and a handy intrinsics guide for people new to SIMDS I
  found online [http://const.me/articles/simd/simd.pdf](
  http://const.me/articles/simd/simd.pdf).
