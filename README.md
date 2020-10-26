# RNA-seq analysis of SMG5-SMG7 in NMD
___
Code and scripts for the RNA-seq analysis of the project: __Nonsense-mediated mRNA decay relies on "two-factor authentication" by SMG5-SMG7__

## Graphical abstract

<img src="https://github.com/boehmv/SMG5-SMG7/blob/main/2FA.png?raw=true" max-height="300">

## Features / Requirements
* Complete analysis of RNA-seq data (ArrayExpress: [E-MTAB-9330](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9330/); provided in FASTQ format), mapped to Gencode v33 / GRCh38.primary_assembly supplemented with SIRVomeERCCome (from Lexogen) using STAR, followed by transcript quantification using Salmon in mapping-based mode with a decoy-aware transcriptome index, finished with analyses of differential gene expression (DGE) via DESeq2, differential transcript usage (DTU) via IsoformSwitchAnalyzeR, alternative splicing (AS) via LeafCutter and intron retention (IR) via IRFinder.
* The main Bash script [CRSA_V006.sh](https://github.com/boehmv/SMG5-SMG7/blob/main/Code/CRSA_V006.sh) requires a design file specifying reference type, sequencing design (single- ,paired-end), study name, folder locations (srvdir for raw file locations, mydir for analyses output) and location of the experiment file. Please see the provided design file example for more information [design.txt](https://github.com/boehmv/SMG5-SMG7/blob/main/Code/design.txt). The tab-delimited experiment file specifies sample name and condition [experiment.txt](https://github.com/boehmv/SMG5-SMG7/blob/main/Code/experiment.txt). Please see the comments in CRSA_V006.sh for further help or run CRSA_V006.sh -h 
* To run/reproduce the complete analysis script, please make sure you have the following tools installed and configured if required:
  * [STAR](https://github.com/alexdobin/STAR) - version 2.7.3a was used for the analyses - with genome indices generated using GRCh38.primary.SIRVomeERCCome.fa and gencode.v33.SIRVomeERCCome.annotation.gtf (both reference files can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). The following code was used for genome index generation: 
  ```
  STAR   --runMode genomeGenerate   --runThreadN 15   --genomeDir /home/volker/reference/gencode.v33.SIRVomeERCCome   --genomeFastaFiles /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa      --sjdbGTFfile /home/volker/reference/Gencode/gencode.v33.SIRVomeERCCome.annotation.gtf   --sjdbOverhang 99
  ```
  * [Alfred](https://github.com/tobiasrausch/alfred) - version v0.2.1 was used for the analyses
  * [samtools](http://www.htslib.org/) - version 1.9 (using htslib 1.9) was used for the analyses
  * [IGV tools](http://software.broadinstitute.org/software/igv/download) - version 2.8.0 was used for the analyses - make sure you have the gencode.v33.SIRVomeERCCome.chrom.sizes file (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)) located in /PATH/TO/IGV_2.8.0/lib/genomes
  * [Salmon](https://github.com/COMBINE-lab/salmon) - version v1.3.0 was used for the analyses - with an index generated using gentrome.v33.SIRV.ERCC.fa.gz and decoys.txt (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). A separate conda environment was created for Salmon. The following code was used for index generation: 
  ```
  salmon index -t /home/volker/reference/Gencode/gentrome.v33.SIRV.ERCC.fa.gz -d /home/volker/reference/Gencode/decoys.txt -p 12 -i /home/volker/reference/Transcriptome/gencode.v33.SIRVomeERCCome --gencode
  ```
  * [DESeq2](https://github.com/mikelove/DESeq2) - version 1.28.1 was used for the analyses - please see [R_sessions](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/R_sessions) for details on R version and other installed packages. The tx2gene file used for the analyses can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)
  * [IsoformSwitchAnalyzeR](https://github.com/kvittingseerup/IsoformSwitchAnalyzeR) - version 1.10.0 was used for the analyses - please see [R_sessions](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/R_sessions) for details on R version and other installed packages.
  * [LeafCutter](https://github.com/davidaknowles/leafcutter) - version v0.2.7 was used for the analyses - please see [R_sessions](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/R_sessions) for details on R version and other installed packages. **NOTE:** small changes in the /scripts of LeafCutter maintained gene IDs from Gencode (changed in [gtf_to_exons_vb.R](https://github.com/boehmv/SMG5-SMG7/blob/main/Code/Tools/leafcutter/scripts/gtf_to_exons_vb.R) and [leafcutter_ds.R](https://github.com/boehmv/SMG5-SMG7/blob/main/Code/Tools/leafcutter/scripts/leafcutter_ds.R)
  * [IRFinder](https://github.com/williamritchie/IRFinder) - version 1.2.6 was used for the analyses - **NOTE:** IRFinder results were not used in the publication
  * [FastQC](https://github.com/s-andrews/FastQC) - version 0.11.9 was used for the analyses
  * [MultiQC](https://github.com/ewels/MultiQC) - version v1.8  was used for the analyses
* All analyses were performed on a 16-core (2x Intel(R) Xeon(R) CPU E5-2687W v2 @ 3.40GHz) workstation with 128 GB RAM running Ubuntu 18.04.4 LTS
* Please make sure to change installation and file paths in the respective scripts

## Log files
Please see [here](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/Log_files) to access the log files from the complete analysis (run between 04.08.2020 - 10.08.2020)

## Quality control
Please see [here](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/QC) to access the MultiQC result file for all samples

## R_sessions
Please see [here](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/R_sessions) to access the session info for the individual R scripts

## Individual scripts
The specialized scripts called by the main CRSA_V006.sh script can be found [here](https://github.com/boehmv/SMG5-SMG7/tree/main/Code/Tools). 

## Feedback / Questions
Feedback is welcome! For any question, please email: boehmv@uni.koeln.de or [create an issue](https://github.com/boehmv/SMG5-SMG7/issues)

## Citation
Volker Boehm, Sabrina Kueckelmann, Jennifer V. Gerbracht, Thiago Britto-Borges, Janine Altmüller, Christoph Dieterich and Niels H. Gehring (2020) __Nonsense-mediated mRNA decay relies on “two-factor authentication” by SMG5-SMG7__. 
bioRxiv 2020.07.07.191437; doi: https://doi.org/10.1101/2020.07.07.191437
