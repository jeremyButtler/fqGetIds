#!/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TOC:
#    sec-1: Variable declerationS, check input, check output
#    sec-2: Get read counts
#    sec-3: Loop though each percentage (hard coded) & build filter file
#    sec-4: Run blank cases to ensure to create more stability
#    sec-5: Run time trials
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerationS, check input, check output
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

fastqStr="$1";
numRepInt="$2";    # number of replicates to do
illBl=$3;           # 1 for Illumina; else nanopore
readsInStartFqI=0; # Number of reads in unsampled fastq
readsInFastqI=0; # holds number of reads in fastqFile
intRep=1;
numReadsInt=0;   # number reads extracted
maxThreads=7;    # max number of threads to test
statsFileStr="fqGetIds-seqkit-benchmark.tsv";
#sedCmdStr="p;n;n;n;n;n;"; # for Illumina (fqGetIds needs)
sedCmdStr="p;n;n;n;"; # for Nanopore
scriptDir="$(dirname $0)";

if [[ ! "$illBl" -eq 0 ]]; then
  sedCmdStr="p;n;n;n;n;n;";
fi # If doing Illumina; then fqGetIds has 6 line output

if [[ ! -f "$fastqStr" ]]; then
    printf \
        "No fastq provided (first argument)\n";
fi # if no ids provided

if [[ "$numRepInt" == "" ]]; then
    numRepInt=10;
fi # if user did not provide a number of replicates to do

if [[ ! -f "$statsFileStr" ]]; then
    { # merge printf commands
        printf "Program\tfastqFile\treadsInFastq\textractionSize";
        printf "\treadsExtracted\tReplicate\tThreadsOrVect\tElapsedTime";
        printf "\tCPUKernTime\tCPUserTime\tMaxResidentMemoryInKb";
        printf "\tPercentCPU\tfileStatus\n";
    } > "$statsFileStr";
fi # if need to create the stats file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Get read counts
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

readsInStartFqI="$(
    awk '
        { # MAIN
            intReads++;         # header
            numSeqLinesInt = 0; # number of lines in sequence entry

            # Nanopore entries are one line, but Illumina likes to use two line
            # So I need a more flexible system
            while($0 !~ /^\+/)
            { # while are are sequence entries
                numSeqLinesInt++;
                getline;
            } # while are are sequence entries

            getline;            # move to fastq entry

            for(intQ = 0; intQ < numSeqLinesInt; intQ++)
                getline; # move to header entry
        }; # MAIN block

        END{print intReads;};
    ' < "$fastqStr" \
 )"; # get the number of reads in the fastq file

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3:
#   - Hot test
#   Loop though each percentage (hard coded) & build filter file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [[ -f "$scriptDir/benchmark-offcache.filt" ]]; then
    rm "$scriptDir/benchmark-offcache.filt"; # make sure awk has new temporary filke
fi # make sure not test file

for intPer in 1 10 100 1000 2250 3500 5000 7500 10000; do
# For the tested percentages

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Sec-5: Run time for se kit trials
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    intRep=1;


    while [[ "$intRep" -le "$numRepInt" ]]; do
    # For all desiered replicates
        iThread=1;

        while [[ "$iThread" -le "$maxThreads" ]]; do
        # while have threads to test
            # Extract 75% of reads for each test
            readsInFastqI="$(
                awk -f "$scriptDir/sampleReads.awk" \
                    -v numReadsToKeepI="$(((75 * readsInStartFqI) /100))" \
                    -v readsInFqI="$readsInStartFqI"\
                    -v offsetI=0 \
                    -v prefixStr="$scriptDir/benchmark-offcache" \
                  < "$fastqStr" \
            )"; # make filter file & get counts

            numReadsInt="$(
                awk -f "$scriptDir/sampleReadIds.awk" \
                    -v numIdsToKeepI="$(((intPer * readsInFastqI) /10000))" \
                    -v readsInFqI="$readsInFastqI" \
                    -v offsetI=0 \
                    -v prefixStr="$scriptDir/benchmark-offcache" \
                  < "$scriptDir/benchmark-offcache.fastq" \
            )"; # make filter file & get counts

            head \
                -n "$(((75 * readsInStartFqI) /25))" \
                < "$fastqStr" \
                > "$scriptDir/delete-file.fastq";
            rm "$scriptDir/delete-file.fastq";

            tail \
                -n "$(((75 * readsInStartFqI) /25))" \
                < "$fastqStr" \
                > "$scriptDir/delete-file.fastq";
            rm "$scriptDir/delete-file.fastq";

            head \
                -n "$(((75 * readsInStartFqI) /25))" \
                < "$fastqStr" \
                > "$scriptDir/delete-file.fastq";
            rm "$scriptDir/delete-file.fastq";

           # Set up meta data for the stats file
           metaStr="$fastqStr	$readsInFastqI	$numReadsInt";

           /usr/bin/time \
               -f "%e\t%S\t%U\t%M\t%P" \
               -o "$scriptDir/tmp-time.tsv" \
             seqkit \
                grep \
                -f "$scriptDir/benchmark-offcache.filt" \
                -j "$iThread" \
                "$scriptDir/benchmark-offcache.fastq" \
             > "$scriptDir/tmp-test-benchmark.fastq"

           timeStr="$(cat "$scriptDir/tmp-time.tsv")";

           numReadsI="$(
               sed -n 'p;n;n;n;' "$scriptDir/tmp-test-benchmark.fastq" | wc -l
           )"; # find number of reads extracted
               # seqkit converts the Illumina fastq 6 read entries into
               # 4 read entries, so do not need a sed command

           printf "seqkit\t%s\t%s\t%s\t%s\t%s\toffCache\n" \
               "$metaStr" \
               "$numReadsI" \
               "$intRep" \
               "$iThread" \
               "$timeStr" \
             >> "$statsFileStr";
           
           iThread=$((iThread + 1));

           rm "$scriptDir/benchmark-offcache.fastq";
           rm "$scriptDir/benchmark-offcache.filt";
        done # while have threads to test

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-6: Get time and stats for fqGetIds run with AVX complie
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        # Extract 75% of reads for each test
        readsInFastqI="$(
            awk -f "$scriptDir/sampleReads.awk" \
                -v numReadsToKeepI="$(((75 * readsInStartFqI) /100))" \
                -v readsInFqI="$readsInStartFqI"\
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$fastqStr" \
        )"; # make filter file & get counts

        numReadsInt="$(
            awk -f "$scriptDir/sampleReadIds.awk" \
                -v numIdsToKeepI="$(((intPer * readsInFastqI) /10000))" \
                -v readsInFqI="$readsInFastqI" \
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$scriptDir/benchmark-offcache.fastq" \
        )"; # make filter file & get counts

        # Run a few blanks to clear the cache
        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        tail \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        # Set up meta data for the stats file
        metaStr="$fastqStr	$readsInFastqI	$numReadsInt";

        /usr/bin/time \
            -f "%e\t%S\t%U\t%M\t%P" \
            -o "$scriptDir/tmp-time.tsv" \
          "$scriptDir/../"fqGetIdsAVX2 \
            -f "$scriptDir/benchmark-offcache.filt" \
            -fastq "$scriptDir/benchmark-offcache.fastq" \
          > "$scriptDir/tmp-test-benchmark.fastq";

        timeStr="$(cat "$scriptDir/tmp-time.tsv")";
        numReadsI="$(
            sed -n "$sedCmdStr" "$scriptDir/tmp-test-benchmark.fastq" | wc -l
        )"; # find number of reads extracted

        printf "fqGetIds\t%s\t%s\t%s\t%s\t%s\toffCache\n" \
            "$metaStr" \
            "$numReadsI" \
            "$intRep" \
            "AVX2" \
            "$timeStr" \
          >> "$statsFileStr";
       
        rm "$scriptDir/benchmark-offcache.fastq";
        rm "$scriptDir/benchmark-offcache.filt";

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-7: Get time and stats for fqGetIds with SSE compile
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        # Extract 75% of reads for each test
        readsInFastqI="$(
            awk -f "$scriptDir/sampleReads.awk" \
                -v numReadsToKeepI="$(((75 * readsInStartFqI) /100))" \
                -v readsInFqI="$readsInStartFqI"\
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$fastqStr" \
        )"; # make filter file & get counts

        numReadsInt="$(
            awk -f "$scriptDir/sampleReadIds.awk" \
                -v numIdsToKeepI="$(((intPer * readsInFastqI) /10000))" \
                -v readsInFqI="$readsInFastqI" \
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$scriptDir/benchmark-offcache.fastq" \
        )"; # make filter file & get counts

        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        tail \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        # Set up meta data for the stats file
        metaStr="$fastqStr	$readsInFastqI	$numReadsInt";

        /usr/bin/time \
            -f "%e\t%S\t%U\t%M\t%P" \
            -o "$scriptDir/tmp-time.tsv" \
          "$scriptDir/../"fqGetIdsSSE \
              -f "$scriptDir/benchmark-offcache.filt" \
              -fastq "$scriptDir/benchmark-offcache.fastq" \
          > "$scriptDir/tmp-test-benchmark.fastq";

        timeStr="$(cat "$scriptDir/tmp-time.tsv")";
        numReadsI="$(
            sed -n "$sedCmdStr" "$scriptDir/tmp-test-benchmark.fastq" | wc -l
        )"; # find number of reads extracted

        printf "fqGetIds\t%s\t%s\t%s\t%s\t%s\toffCache\n" \
            "$metaStr" \
            "$numReadsI" \
            "$intRep" \
            "SSE" \
            "$timeStr" \
          >> "$statsFileStr";

        rm "$scriptDir/benchmark-offcache.fastq";
        rm "$scriptDir/benchmark-offcache.filt";

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-8: Get time and stats for fqGetIds with scalar complie
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        # Extract 75% of reads for each test
        readsInFastqI="$(
            awk -f "$scriptDir/sampleReads.awk" \
                -v numReadsToKeepI="$(((75 * readsInStartFqI) /100))" \
                -v readsInFqI="$readsInStartFqI"\
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$fastqStr" \
        )"; # make filter file & get counts

        numReadsInt="$(
            awk -f "$scriptDir/sampleReadIds.awk" \
                -v numIdsToKeepI="$(((intPer * readsInFastqI) /10000))" \
                -v readsInFqI="$readsInFastqI" \
                -v offsetI=0 \
                -v prefixStr="$scriptDir/benchmark-offcache" \
              < "$scriptDir/benchmark-offcache.fastq" \
        )"; # make filter file & get counts

        # Run a few blanks to clear the cache
        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        tail \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        head \
            -n "$(((75 * readsInStartFqI) /25))" \
            < "$fastqStr" \
            > "$scriptDir/delete-file.fastq";
        rm "$scriptDir/delete-file.fastq";

        metaStr="$fastqStr	$readsInFastqI	$numReadsInt";

        /usr/bin/time \
            -f "%e\t%S\t%U\t%M\t%P" \
            -o "$scriptDir/tmp-time.tsv" \
          "$scriptDir/../"fqGetIdsScalar \
              -f "$scriptDir/benchmark-offcache.filt" \
              -fastq "$scriptDir/benchmark-offcache.fastq" \
          > "$scriptDir/tmp-test-benchmark.fastq";

        timeStr="$(cat "$scriptDir/tmp-time.tsv")";
        numReadsI="$(
            sed -n "$sedCmdStr" "$scriptDir/tmp-test-benchmark.fastq" | wc -l
        )"; # find number of reads extracted

        printf "fqGetIds\t%s\t%s\t%s\t%s\t%s\toffCache\n" \
            "$metaStr" \
            "$numReadsI" \
            "$intRep" \
            "Scalar" \
            "$timeStr" \
          >> "$statsFileStr";

        rm "$scriptDir/benchmark-offcache.fastq";
        rm "$scriptDir/benchmark-offcache.filt";

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Sec-9: Clean up
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        intRep=$((intRep + 1));
    done
done # For the tested percentages

rm "$scriptDir/tmp-time.tsv";

exit;
