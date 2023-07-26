
# Input:
#   -v numReadsToKeepI=1000
#   -v readsInFqI=60000000
#   -v offsetI=0
#   -v prefixStr=out

BEGIN{
  if(offsetI == "") offsetI = 0;
  if(prefixStr == "") prefixStr = "out";

  totalReadsI = 0; # Will hold number of read ids I read in
  modValI = int(readsInFqI / numReadsToKeepI);
}; # BEGIN

{ # MAIN
  if(keptReadsI > numReadsToKeepI) exit;

  totalReadsI++;  # count number of reads I have read in

  if(totalReadsI % modValI == offsetI)
  { # if is a read I want to keep
    # clean up read id for seqkit (fastqGrep can handle this)
    outFILE = prefixStr ".fastq";
    numSeqLinesInt = 0;

    # print out the header
    print $0 >> outFILE;
    getline;

    while($1 !~ /^\+/)
    { # while are are sequence entries
      print $0 >> outFILE;
      numSeqLinesInt++;
      getline;
    } # while are are sequence entries

    # print out the spacer
    print $0 >> outFILE;
    getline;

    # prit out first qscore entry
    print $0 >> outFILE;
    numSeqLinesInt--;

    # print out the Q-score entry
    while(numSeqLinesInt > 0)
    { # For the q-score entry
      getline;
      print $0 >> outFILE;
      numSeqLinesInt--;
    } # For the q-score entry

    keptReadsI++;  # keep track of how many ideas I kept
    next;
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

END{print keptReadsI}; # make sure user gets counts
