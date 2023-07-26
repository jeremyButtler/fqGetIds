#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SOP: Start of program
#   o sec-0: Includes and functions
#   o sec-1: Variable declerations
#   o sec-2: Do some clean up of the dataset
#   o sec-3: Graph the elapsed time for Nanopore
#   o sec-4: Graph the elapsed time for Illumina
#   o sec-5: Graph the memory use (maximum resident)
#   o sec-6: Graph percent cpu used
#   o sec-7: Graph the pi Illumina elapsed time separately
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-0: Indlues and functions
#   o sec-0 sub-1: Included libraries
#   o sec-0 sub-2: saveGraph function
#     - Wrapper for ggsave, here to allow easy changes
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-0 Sub-1: Included libraries
#**********************************************************

library("ggplot2")        # graphing
library("data.table")     # apply functions to dataframe gourps
library("ggpubr")         # theme
library("viridis")        # color pallete

#**********************************************************
# Sec-0 Sub-2:
#  Output: Saves output graph as an svg
#**********************************************************
saveGraph = function(
  nameStr           # name of graph to save
) # Use: saves a graph using ggsave function
{ # saveGraph
    ggsave(paste(nameStr, ".svg", sep = ""), 
           device = "svg", # save as tiff file
           dpi = 300,
    ); # ggsave (save the graph)
} # saveGraph

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-1: Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

fileStr="fqGetIds-seqkit-benchmark.tsv"
benchDT = setDT(read.csv(fileStr, header = TRUE, sep = "\t"));
    # setDT converts to data table (easier to manipulate)
graphDT = NULL; # so can manipulate benchDT safely

graphObj = NULL;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-2: Do some clean up of the dataset
#   o sec-2 sub-1:
#     - remove fqGetIdsMem (suspect wrong compile)
#   o sec-2 sub-2:
#     - Number conversion to get % reads & memory usage in
#       Mb
#   o sec-2 sub-3:
#     - Rename catagories to have clearer names
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-2 Sub-1:
#  - remove fqGetIdsMem (suspect wrong compile)
#**********************************************************

benchDT = benchDT[benchDT$Program != "fqGetIdsMem",];
benchDT$test =
  paste(benchDT$Program, benchDT$Threads, sep = "-");

#**********************************************************
# Sec-2 Sub-2:
#  - Number conversions to get % reads & memory usage in Mb
#**********************************************************

# Get the percentage of reads extracted from the fastq
benchDT$percReads =
  100 * (benchDT$extractionSize / benchDT$readsInFastq);

# Get maximum resident memory usage in megabytes
benchDT$maxResMemInMb =
  benchDT$MaxResidentMemoryInKb / 1000;

# Remove the % symbol in PercentCPU column (so is numeric)
benchDT$PercentCPU =
  as.numeric(sub("%", "", benchDT$PercentCPU));
    # as.numeric(): convert string to number
    # sub: remove one % symbol

#**********************************************************
# Sec-2 Sub-3:
#  - Rename catagories to have clearer names
#**********************************************************

# Set up IO limit names
#benchDT$Computer =
#  sub("fastIO", "Fast IO", benchDT$Computer);
#benchDT$Computer =
#  sub("piSlowIO", "IO limited (pi 1B)", benchDT$Computer);
# This was for when I was able to run seqkit on a pi 1B.
#   Saddly the version of seqkit on the pi is 0.15 and
#   needs to be updated (manual install)

benchDT$tech = benchDT$fastqFile;

# Reset the Illumina file names to the sequencer
benchDT$tech =
   sub(
       "../../seqkit-benchmark-data/dataset_C.fq",
       #"Illumina fast=9186045/pi=890841 reads",
       "Illumina = 6889534 reads",
       benchDT$tech
); # This is only the offcache tests

benchDT$tech =
   sub(
       "../../seqkit-benchmark-data/75PercIll.fastq",
       "Illumina = 6889534 reads",
       benchDT$tech
); # This was the on cache tests


benchDT$tech =
   sub(
       "../../AllUP9.fastq",
       "Nanopore = 4523546 reads",
#       "Nanopore 6031394/890841 reads",
       benchDT$tech
); # nanopore on cache tests

benchDT$tech =
   sub(
       "../../75PercAllUP9.fastq",
       "Nanopore = 4523546 reads",
#       "Nanopore 6031394/890841 reads",
       benchDT$tech
); # nanopore on cache tests

# This was for the pi
#benchDT$tech =
#    sub(
#        "UP9_20191212_LSK109_IVM.barcode10.qcreads.fastq",
#        "Nanopore fast=6031394/pi=890841 reads",
#        benchDT$tech
#); # Pi 1B (slow IO) Nanopore data
#
#benchDT$tech =
#    sub(
#        "Ill-890841-reads.fastq",
#        "Illumina fast=9186045/pi=890841 reads",
#        benchDT$tech
#); # Pi 1B (slow IO) Illumina data

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-3:
#  - Graph the elapsed time for Nanopore
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# make the graph for the on cache data
graphObj =
  ggplot(
    data =
      benchDT[
        benchDT$fileStatus != "offCache" &
        benchDT$readsInFastq < 6889534 &
        benchDT$ThreadsOrVect != 2 &
        benchDT$ThreadsOrVect != 4 &
        benchDT$ThreadsOrVect != 5 &
        benchDT$ThreadsOrVect != 6
        ,
       ],
    aes(y = ElapsedTime, x = percReads)
);
            
graphObj =
    graphObj +
    geom_vline(
        aes(xintercept = 50),          # make verical line at 50% mark
        linetype = "dashed"            # make a dashed line
    ) +
    geom_point(
        aes(col = test, shape = test), # color and shape of each point
            # color & shape determined by program & number threads used
        alpha = 0.5,                 # transparency of each point
        size = 4,                    # size of each point
        position = position_jitter(height = 0, width = 3)
             # Add in a jitter (space points out by x-axis (width))
    ) + # Plot the points
    facet_grid(
        cols = vars(tech),  # make columns (by sequencing tech)
        rows = vars(fileStatus),   # rows to separate IO slow from IO fast
        scales = "free_y"        # Rows have different y-axis scales
    ) +
    scale_color_viridis(             # add color
        direction = -1,              # color in oposite direction
        option = "D",                # Use the default color pallet
        discrete = TRUE,             # Is discrete data (not continuous)
        name = "Program-threads"     # Name of legend
    ) +
    scale_shape_discrete(            # Add shapes to my data
        name = "Program-threads"     # Name of legend
    ) +        
    ylab("Time in seconds") +        # y-axis title
    xlab("Percentage of reads extracted") + # x-axis title
    theme_pubr() +                   # theme I like to use
    theme(axis.text.x = element_text(angle = 90)); #rotate x-axis labels
        # needs to come after theme_pubr(), otherwise pubr overwrites

saveGraph("fqGetIds-nano-time");
#dev.off(); # clear the last plot

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-4:
#  - Graph the elapsed time for Illumina
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# make the graph for the on cache data
graphObj =
  ggplot(
    data =
      benchDT[
        benchDT$fileStatus != "offCache" &
        benchDT$readsInFastq == 6889534 &
        benchDT$ThreadsOrVect != 2 &
        benchDT$ThreadsOrVect != 4 &
        benchDT$ThreadsOrVect != 5 &
        benchDT$ThreadsOrVect != 6
        ,
       ],
    aes(y = ElapsedTime, x = percReads)
);
            
graphObj =
    graphObj +
    geom_vline(
        aes(xintercept = 50),          # make verical line at 50% mark
        linetype = "dashed"            # make a dashed line
    ) +
    geom_point(
        aes(col = test, shape = test), # color and shape of each point
            # color & shape determined by program & number threads used
        alpha = 0.5,                 # transparency of each point
        size = 4,                    # size of each point
        position = position_jitter(height = 0, width = 3)
             # Add in a jitter (space points out by x-axis (width))
    ) + # Plot the points
    facet_grid(
        cols = vars(tech),  # make columns (by sequencing tech)
        rows = vars(fileStatus),   # rows to separate IO slow from IO fast
        scales = "free_y"        # Rows have different y-axis scales
    ) +
    scale_color_viridis(             # add color
        direction = -1,              # color in oposite direction
        option = "D",                # Use the default color pallet
        discrete = TRUE,             # Is discrete data (not continuous)
        name = "Program-threads"     # Name of legend
    ) +
    scale_shape_discrete(            # Add shapes to my data
        name = "Program-threads"     # Name of legend
    ) +        
    ylab("Time in seconds") +        # y-axis title
    xlab("Percentage of reads extracted") + # x-axis title
    theme_pubr() +                   # theme I like to use
    theme(axis.text.x = element_text(angle = 90)); #rotate x-axis labels
        # needs to come after theme_pubr(), otherwise pubr overwrites

saveGraph("fqGetIds-Ill-time");
#dev.off(); # clear the last plot

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5:
#  - Graph the memory use (maximum resident)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# make the graph
graphObj=
  ggplot(
    data =
      benchDT[
        benchDT$fileStatus != "offCache" &
        benchDT$ThreadsOrVect != 2 &
        benchDT$ThreadsOrVect != 4 &
        benchDT$ThreadsOrVect != 5 &
        benchDT$ThreadsOrVect != 6
        ,
       ],
    aes(y = maxResMemInMb, x = percReads)
);
            
graphObj =
    graphObj +
    geom_point(
        aes(col = test, shape = test), # color and shape of each point
            # color & shape determined by program & number threads used
        alpha = 0.5,                 # transparency of each point
        size = 4,                    # size of each point
        position = position_jitter(height = 0, width = 3)
             # Add in a jitter (space points out by x-axis (width))
    ) + # Plot the points
    facet_grid(
        cols = vars(tech),  # make columns (by sequencing tech)
        rows = vars(fileStatus),   # rows to separate IO slow from IO fast
        scales = "free_y"        # Rows have different y-axis scales
    ) + # make rows (by computer used)
    scale_color_viridis(             # add color
        direction = -1,              # color in oposite direction
        option = "D",                # Use the default color pallet
        discrete = TRUE,             # Is discrete data (not continuous)
        name = "Program-threads"     # Name of legend
    ) +
    scale_shape_discrete(            # Add shapes to my data
        name = "Program-threads"     # Name of legend
    ) +        
    ylab("Maximum resident memory usage in megabytes") + # y-axis title
    xlab("Percentage of reads extracted") + # x-axis title
    theme_pubr() +                   # theme I like to use
    theme(axis.text.x = element_text(angle = 90)); #rotate x-axis labels
        # needs to come after theme_pubr(), otherwise pubr overwrites

saveGraph("fqGetIds-memory");
#dev.off(); # clear the last plot

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-5: Graph percent cpu used
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# make the graph
graphObj= ggplot(
    data =
      benchDT[
        benchDT$fileStatus != "offCache" &
        benchDT$ThreadsOrVect != 2 &
        benchDT$ThreadsOrVect != 4 &
        benchDT$ThreadsOrVect != 5 &
        benchDT$ThreadsOrVect != 6
        ,
       ],
    aes(y = PercentCPU, x = percReads)
);
            
graphObj =
    graphObj +
    geom_vline(
        aes(xintercept = 50),          # make verical line at 50% mark
        linetype = "dashed"            # make a dashed line
    ) +
    geom_point(
        aes(col = test, shape = test), # color and shape of each point
            # color & shape determined by program & number threads used
        alpha = 0.5,                 # transparency of each point
        size = 4,                    # size of each point
        position = position_jitter(height = 0, width = 3)
             # Add in a jitter (space points out by x-axis (width))
    ) + # Plot the points
    facet_grid(
        cols = vars(tech),  # make columns (by sequencing tech)
        rows = vars(fileStatus),   # rows to separate IO slow from IO fast
        scales = "free_y"        # Rows have different y-axis scales
    ) + # make rows (by computer used)
    scale_color_viridis(             # add color
        direction = -1,              # color in oposite direction
        option = "D",                # Use the default color pallet
        discrete = TRUE,             # Is discrete data (not continuous)
        name = "Program-threads"     # Name of legend
    ) +
    scale_shape_discrete(            # Add shapes to my data
        name = "Program-threads"     # Name of legend
    ) +        
    ylab("Percent of CPU used") + # y-axis title
    xlab("Percentage of reads extracted") + # x-axis title
    ylim(0, NA) +                 # change y-axis limits
        # 0, NA: set lowest limit to 0, & do not change max limit (NA)
    theme_pubr() +                   # theme I like to use
    theme(axis.text.x = element_text(angle = 90)); #rotate x-axis labels
        # needs to come after theme_pubr(), otherwise pubr overwrites

saveGraph("fqGetIds-cpu");
#dev.off(); # clear the last plot

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Sec-7:
##  - Graph the pi Illumina elapsed time separately
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
## make the graph
#graphObj =
#  ggplot(
#    data = benchDT[benchDT$fastqFile == "Ill-890841-reads.fastq",],
#    aes(y = ElapsedTime, x = percReads)
#); # graph data
#            
#graphObj =
#    graphObj +
#    geom_vline(
#        aes(xintercept = 50),          # make verical line at 50% mark
#        linetype = "dashed"            # make a dashed line
#    ) +
#    geom_point(
#        aes(col = test, shape = test), # color and shape of each point
#            # color & shape determined by program & number threads used
#        alpha = 0.5,                 # transparency of each point
#        size = 4,                    # size of each point
#        position = position_jitter(height = 0, width = 3)
#             # Add in a jitter (space points out by x-axis (width))
#    ) + # Plot the points
#    scale_color_viridis(             # add color
#        direction = -1,              # color in oposite direction
#        option = "D",                # Use the default color pallet
#        discrete = TRUE,             # Is discrete data (not continuous)
#        name = "Program-threads"     # Name of legend
#    ) +
#    scale_shape_discrete(            # Add shapes to my data
#        name = "Program-threads"     # Name of legend
#    ) +        
#    ylab("Time in seconds") +        # y-axis title
#    xlab(
#      "Pi 1B Illumina, percentage of reads extracted from 890841 reads"
#    ) + # x-axis title
#    theme_pubr() +                   # theme I like to use
#    theme(axis.text.x = element_text(angle = 90)); #rotate x-axis labels
#        # needs to come after theme_pubr(), otherwise pubr overwrites
#
#saveGraph("fqGetIds-seqkit-full-bench--1b-Ill-elapsed-time");
##dev.off(); # clear the last plot

