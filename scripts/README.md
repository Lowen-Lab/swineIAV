# Scripts necessary for reproducing the analyses from our manuscript

To reproduce our analyses, the scripts in this directory should be used in the following order:
1. After collecting the data from the SRA under bioprojectID PRJNA1051292, place those files in a folder called "raw_data". Then retrieve the reference sequence files in the "manuscript_data" folder in this github repository and place them in a folder called "refs"

2. run "map_reads.py". This script takes the raw reads as input and uses bbtools and other software to screen the reads for adapters and split reads according to IAV segment before mapping those reads to references sequences.
   * It is written to run all samples in parallel using all of the cores you make available to it.
   * Analysis can be interrupted at anytime and it will resume without repeating any steps. It has safeguards to prevent incompletely processed samples from being included in subsequent analyses by confirming that the exit code from all steps in the pipeline return no error codes, and by only transferring files from a temporary working directory to their permanent directory after all analysis from that stage in the pipeline is complete.

3. Once all reads have been processed and mapped, run "SAM_parse.py". This script parses the .sam files generated in the previous step to collate properties of all alleles at every site detected in the alignment.
   * This script is written similarly to "map_read.py" in that is processes all samples in parallel and performs repeated sanity checks to confirm analysis was complete prior to reporting final values
   * The output of this script is a ".substitutions.txt" file that contains summary statistics on all nucleotides at all sites in each segment, including:
      * nucleotide frequency
      * average base quality
      * average mapping quality
      * average read length of reads containing each nucleotide at that site
      * average distance of that allele to the nearest end of the read (erroneous iSNVs are often predominantly found in the first (or last) 5-10 bases of illumina reads despite having very high phred scores)
      * average number of mismatches to the reference sequences that reads containing that allele have
      * average number of indels to the reference sequences that reads containing that allele have

4. Once the SAM_parse.py script has finished processing all samples, you can now run "SNV_classify.py". This script generates a large number of tab delimited tables that were directly used in all analyses and figures from our manuscriptvalues
   * This script became very purpose-built to the swine IAV manuscript, such as by performing checks to make sure samples are linked to valid pigIDs, manually curated genotype entries, and other things added to handle the exceptions, structure, and quirks unique to our dataset. Accordingly, it will not be able to process similar WGS datasets without considerable revisions. If you would like to process a similar dataset using our workflow, please contact the administrator of this account.
   * The major outputs from this script are:
      * iSNV tables that reflect frequency in pigs over time (or labeled as not detected where time points contained insufficient coverage or sequence/mapping quality to determine if the iSNV was detected or not)
      * Synonymous and non-synonymous iSNV frequency and ratios necessary for comparing observed ratios to a neutral expectation (calculated by this script using the consensus sequences it creates)
      * summary files with coverage, base quality, mapping quality, etc for the major and minor alleles in all samples for emperically setting the thresholds for iSNV detection, and for diagnosing any issues with how samples are being processed or mapped

5. To generate the haplotype network and detect dSNPs using different frequency, depth, and other cutoffs, run "MSA_to_haplotype_network.py".
   * This script requires the python module NetworkX, which can be easily installed with conda.
   * All of the tools for generating haplotype networks that we tried did not handle input data with Ns from low coverage regions, so we wrote this script to both identify the lineage defining SNPs used in our genotyping method, but also generate the haplotype network figure
   * The graphical ouput of this script is a .graphml file containing the minimum spanning tree networks of each samples' segments based on pairwise sequence distance (the haplotype network), which can be opened in software like gephi to generate publication quality figures.
