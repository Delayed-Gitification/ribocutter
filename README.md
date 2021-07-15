# Ribocutter

Ribocutter is a tool for designing sgRNA templates to deplete unwanted abundant contaminants from a Ribo-seq (or similar) high-throughput sequencing library. Using an input fastq (that MUST be adaptor trimmed first!), it automatically designs ready-to-order sgRNA templates for use with the NEB single-tube sgRNA synthesis kit. See our paper for more information on the actual experimental protocol

# Usage:

Ribocutter has two required arguments: an input fastq.gz file, and the name of the output. Additionally, you can specificy the exact number of guides you'd like to order and, if you have particularly short contaminating sequences, you can add the flanking 3' and 5' sequences too (which will ideally contain a PAM motif). Note that adding these flanking sequences will increase the risk of off-target depletion - keep them short!

Ribocutter will determine the fraction of the sequencing library that is targeted by your guides. As a rule of thumb, if f is the fraction of reads targeted by guides, the fold-increase in fraction of useful reads will be 1/(1-f). So, if f = 0.67 (67%), then the increase will be 1/(1-0.67) = 3-fold increase. If your guides only target a small fraction of the library (significantly less than 50%) then it's probably not worth your time to apply the protocol as the fold increase will be fairly small - of course you can increase the fraction of reads targeted by increasing the number of guides used (oligos are cheap anyway...)

You can also include a fasta of background sequences that you do not wish to target. However, while this seems like a good option, I would actually advise against this - low-level off-target depletion is a small price to pay for a large increase in useful reads, and your libraries are treated after multiplexing anyway, so this depletion should not lead to erroneous detection of differentially translated genes.

Additionally, you can use ribocutter to produce useful statistics about copy numbers of abundant sequences in your library, by using the --save_stats option.

# Installation:

The easiest way to install is to use pypi (pip install ribocutter). Ribocutter is also on Bioconda (though for latest releases check pip!). If you are so inclined you could also download the python script and run it directly.
