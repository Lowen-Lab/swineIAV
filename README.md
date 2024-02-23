# Swine IAV whole genome sequence analysis
Welcome to the repository for the supporting information and files necessary to reproduce the analyses contained in the paper VanInsberghe et al (2024) "Genetic drift and purifying selection shape within-host influenza A virus populations during natural swine infections"

## How this repository is organized

### Manuscript data
This directory contains all of the files that you will need to download to reproduce the analyses and generate the figures presented in our publication. A more detailed description of those files is contained within.
   * The only data that you will need that is not in this repository is the sequence data, located in the NCBI SRA under the bioproject ID PRJNA1051292.
### Scripts
This directory contains all of the custom python scripts generated for this publication. 
   * We performed all analyses on a linux computing cluster, but it can theoretically be run on any Windows, Mac, or Linux computer that the following dependencies can be installed on:
      1. Python3
         * We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to install python and the required packages
      2. [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
         * Note that the bbtools available on anaconda is not currently a functional copy, so it needs to be downloaded from the JGI website linked above
      3. Java (to run the BBTools)
      4. Bash (also required to run BBTools)
         * This is available by default on Mac and Linus operating systems
         * For windows users, we recommend using [GitBash](https://git-scm.com/downloads). In our testing, WSL was not a reliable method to perform the analysis.
           * To use GitBash, the user will need to edit their *".bashrc" file in *"C:/Users/username"* to GitBash to add the location of their conda profile, the location of their python3 installation, and add the location of their bbtools installation to the GitBash path
              * An example of how you might edit your .bashrc file:
                 . /c/ProgramData/miniconda3/etc/profile.d/conda.sh
                 alias python="winpty /c/ProgramData/miniconda3/python.exe"
                 PATH=$PATH:"/c/ProgramData/bbmap"
      4. Python modules
         * numpy
         * matplotlib
         * multiprocessing
         * joblib
         * networkx
