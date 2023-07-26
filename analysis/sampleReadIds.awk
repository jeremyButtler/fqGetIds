
# Input:
#   -v numIdsToKeepI=1000
#   -v readsInFqI=60000000
#   -v offsetI=0
#   -v prefixStr=out

BEGIN{
  if(offsetI == "") offsetI = 0;
  if(prefixStr == "") prefixStr = "out";

  totalIdsI = 0; # Will hold number of read ids I read in
  modValI = int(readsInFqI / numIdsToKeepI);
}; # BEGIN

{ # MAIN
  if(keptIdsInt > numIdsToKeepI) exit;

  totalIdsI++;  # count number of reads I have read in

  if(totalIdsI % modValI == offsetI)
  { # if is a read I want to keep
    # clean up read id for seqkit (fastqGrep can handle this)
    sub(/^@/, "", $0); # remove @ marking fastq header
    sub(/ .*/, "", $0); # remove extra information in header

    keptIdsInt++;  # keep track of how many ideas I kept
    outFILE = prefixStr ".filt";

    print $0 >> outFILE;  # print out the header
  } # if is a read I want to keep

  numSeqLinesInt = 0;
  getline;            # get off the header

  # Nanopore entries are one line, but Illumina has used
  # two lines for entries. So I need a more flexible system
  while($0 !~ /^\+/)
  { # while are are sequence entries
    numSeqLinesInt++;
    getline;
  } # while are are sequence entries

  # move to header entry
  for(intQ = 0; intQ < numSeqLinesInt; intQ++) getline;
}; # MAIN block

END{print keptIdsInt}; # make sure user gets counts
