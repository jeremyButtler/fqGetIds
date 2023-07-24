CC=cc
PREFIX=/usr/local/bin
NAME=fqGetIds
CFLAGS=\
  -static \
  --std=c99 \
  -O3

FILES=\
  vectorWrap.c \
  fqGetIdsFqFun.c \
  fqGetIdsStructs.c \
  fqGetIdsAVLTree.c \
  fqGetIdsHash.c \
  fqGetIdsSearchFq.c \
  fqGetIds.c

# Program to extract fastq entries by a set of read ids
all:
	$(CC) $(CFLAGS) $(FILES) -o $(NAME) || \
   gcc   $(CFLAGS) $(FILES) -o $(NAME) || \
   egcc  $(CFLAGS) $(FILES) -o $(NAME) || \
   cc    $(CFLAGS) $(FILES) -o $(NAME)

sse:
	$(CC) $(CFLAGS) -DSSE $(FILES) -o $(NAME) || \
   gcc   $(CFLAGS) -DSSE $(FILES) -o $(NAME) || \
   egcc  $(CFLAGS) -DSSE $(FILES) -o $(NAME) || \
   cc    $(CFLAGS) -DSSE $(FILES) -o $(NAME)

avx2:
	$(CC) $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o $(NAME) || \
   gcc   $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o $(NAME) || \
   egcc  $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o $(NAME) || \
   cc    $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o $(NAME)

avx512:
	$(CC) $(CFLAGS) -DAVX512 -mavx512bw $(FILES) -o $(NAME) || \
   gcc   $(CFLAGS) -DAVX512 -mavx512bw $(FILES) -o $(NAME) || \
   egcc  $(CFLAGS) -DAVX512 -mavx512bw $(FILES) -o $(NAME) || \
   cc    $(CFLAGS) -DAVX512 -mavx512bw $(FILES) -o $(NAME)

neon:
	$(CC) $(CFLAGS) -DNEON -mfpu=neon $(FILES) -o $(NAME) || \
   gcc   $(CFLAGS) -DNEON -mfpu=neon $(FILES) -o $(NAME) || \
   egcc  $(CFLAGS) -DNEON -mfpu=neon $(FILES) -o $(NAME) || \
   cc    $(CFLAGS) -DNEON -mfpu=neon $(FILES) -o $(NAME)

install:
	mv $(NAME) $(PREFIX)
	chmod a+x $(PREFIX)/$(NAME)

# For cleaning benchmarking programs
clean:
	rm fqGetIdsScalar fqGetIdsSSE fqGetIdsAVX2 fqGetIdsDebug;

benchmark:
	$(CC) $(CFLAGS) $(FILES) -o fqGetIdsScalar || \
   gcc   $(CFLAGS) $(FILES) -o fqGetIdsScalar || \
   egcc  $(CFLAGS) $(FILES) -o fqGetIdsScalar || \
   cc    $(CFLAGS) $(FILES) -o fqGetIdsScalar
	$(CC) $(CFLAGS) -DSSE $(FILES) -o fqGetIdsSSE || \
   gcc   $(CFLAGS) -DSSE $(FILES) -o fqGetIdsSSE || \
   egcc  $(CFLAGS) -DSSE $(FILES) -o fqGetIdsSSE || \
   cc    $(CFLAGS) -DSSE $(FILES) -o fqGetIdsSSE
	$(CC) $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o fqGetIdsAVX2 || \
   gcc   $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o fqGetIdsAVX2 || \
   egcc  $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o fqGetIdsAVX2 || \
   cc    $(CFLAGS) -DAVX2 -mavx -mavx2 $(FILES) -o fqGetIdsAVX2

# Was for debuging code
egcc:
	egcc -O0 -ggdb -Wall $(FILES) -o $(NAME)Debug
#	#egcc -O3 -ggdb -Wall $(FILES) -o $(NAME)Debug
#	#egcc -O3 -ggdb -Wall -DSSE $(FILES) -o $(NAME)Debug
#	#time ./fqGetIdsDebug -f del.filt -fastq del.fq > del-results.fq
#	egcc -O0 -ggdb -Wall -DSSE $(FILES) -o $(NAME)Debug
#	egdb  \
#     -ex "set autoconfirm off" \
#     -ex "b fqGetIdsSearchFq.c:352" \
#     -ex "b fqGetIdsSearchFq.c:367" \
#     -ex "r -f del.filt -fastq del.fq > del-results.fq" \
#     fqGetIdsDebug
#	rm fqGetIdsDebug
