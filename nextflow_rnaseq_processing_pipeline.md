# nf-core/rnaseq

**Nextflow RNA-seq Processing Pipeline**

*Guide written by Vini Karumuru*

*Last updated: 01/17/2026*

---

# Overview

## Documentation

[rnaseq: Introduction](https://nf-co.re/rnaseq/3.22.2)

- This pipeline performs the following general steps (**optional*):
    1. Quality statistics (FastQC)
    2. UMI extraction (UMI-tools)*
    3. Adapter & quality trimming (Trim Galore!)
    4. Alignment to a reference genome (STAR)
    5. Sorting & indexing of aligned reads (SAMtools)
    6. UMI-based deduplication (UMI-tools)*
    7. Quantification of reads aligning to genes (Salmon)

## CQLS Instructions

[Using Nextflow for scientific computing - Oregon State University Research HPC Documentation](https://docs.hpc.oregonstate.edu/cqls/software/nextflow/)

- Due to memory and set-up requirements, the entirety of this pipeline should be run on the CQLS server
- Nextflow workflows that CQLS has available are stored in `/local/cqls/software/nextflow/assets`
    - In this case, the one to use is `nf-core-rnaseq_3.22.2/`
- Ensure that dot files are updated (only needs to be done once):
    
    [Configuration on the new HPC infrastructure - Oregon State University Research HPC Documentation](https://docs.hpc.oregonstate.edu/cqls/configuration/)
    
    ```bash
    hpcman user update-dotfiles -c all
    exec $SHELL
    ```
    

## Important Considerations

- Some of the parameters selected below for this pipeline are recommended by Lexogen for libraries prepared with their *Quantseq 3’ mRNA-Seq V2 FWD* library preparation kit (with the *UMI Second Strand Synthesis Module*)
    - A different library preparation protocol would require different parameters (see the *Unique Molecular Identifiers* section at [https://nf-co.re/rnaseq/3.21.0/docs/usage/](https://nf-co.re/rnaseq/3.21.0/docs/usage/))
- If a gene count table prior to UMI deduplication is desired, the pipeline must be run twice (once with the UMI deduplication options and once without)

---

# Running the Pipeline

## Set Up

### Creating a STAR Indexed Database

*This step only needs to be completed once for any given reference genome. The STAR index can be re-used for any number of runs.*

1. Download reference genome fasta & GTF files
    - Example: Ensembl 107 dog genome
        - GTF file (from [https://ftp.ensembl.org/pub/release-107/gtf/canis_lupus_familiaris/](https://ftp.ensembl.org/pub/release-107/gtf/canis_lupus_familiaris/)):
            
            ```bash
            wget [https://ftp.ensembl.org/pub/release-107/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.107.gtf.gz](https://ftp.ensembl.org/pub/release-107/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.107.gtf.gz)
            ```
            
        - FASTA file (top level, unmasked genome) (from [https://ftp.ensembl.org/pub/release-107/fasta/canis_lupus_familiaris/dna/](https://ftp.ensembl.org/pub/release-107/fasta/canis_lupus_familiaris/dna/)):
            
            ```bash
            wget https://ftp.ensembl.org/pub/release-107/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz
            ```
            
2. Unzip the genome files
3. Set up and enter a conda environment containing the most up-to-date version of STAR
4. Run the following command via an interactive session or batch job, filling in the last 3 file paths:
    
    ```bash
    STAR 
    --runMode genomeGenerate 
    --runThreadN 128 
    --genomeDir <output_dir> 
    --genomeFastaFiles <path_to_fasta>
    --sjdbGTFfile <path_to_gtf>
    ```
    
5. After running, the output directory should contain the following files:
    
    ```bash
    chrLength.txt
    chrNameLength.txt
    chrName.txt
    chrStart.txt
    exonGeTrInfo.tab
    exonInfo.tab
    geneInfo.tab
    Genome
    genomeParameters.txt
    Log.out
    SA
    SAindex
    sjdbInfo.txt
    sjdbList.fromGTF.out.tab
    sjdbList.out.tab
    transcriptInfo.tab
    ```
    
6. Rezip the genome files

### Samplesheet

1. Create a `.csv` file, that has the following format (column names must be as they appear here):
    
    
    | sample | fastq_1 | fastq_2 | strandedness |
    | --- | --- | --- | --- |
    |  |  |  |  |
2. Fill in the columns with sample information and paths to the input FASTQ files
    - Path to FASTQ files can be absolute paths or must be relative to the location that the run will be started at
    - Leave `fastq_2` empty if the sequencing is single-end
    - Options for `strandedness` are: `unstranded`, `forward`, `reverse`, `auto`
        - If unsure, specify auto for all samples, and the pipeline will detect it for you
    - If a sample was sequenced more than once, it can be listed over multiple rows (repeat sample name) and the pipeline will merge them
    - Do not include empty FASTQ files (i.e. samples with no reads)
3. Example samplesheet:
    
    
    | sample | fastq_1 | fastq_2 | strandedness |
    | --- | --- | --- | --- |
    | Anderson1B | ./input_fastq/Anderson1B_merged.fastq.gz |  | auto |
    | Anderson2B | ./input_fastq/Anderson2B_merged.fastq.gz |  | auto |

### Run Configuration

1. Create a file named exactly `nextflow.config` with the following contents:
    
    ```bash
    singularity.autoMounts = true
    singularity.runOptions = '-B /scratch:/tmp,/scratch'
    process.executor = "slurm"
    process.queue = "biomed"
    process.clusterOptions= '--account=biomed --nodelist=ruby'
    ```
    
2. Create a parameters `.json` file, filling in specified file paths
    - If excluding UMI deduplication:
        
        ```bash
        {
            "input": "<samplesheet.csv>",
            "outdir": "<output_dir>",
            "fasta": "<path_to_fasta>",
            "gtf": "<path_to_gtf>",
            "star_index": "<path_to_star_index>",
            "skip_pseudo_alignment": true,
            "skip_markduplicates": true,
            "skip_bigwig": true,
            "skip_stringtie": true,
            "skip_dupradar": true,
            "skip_qualimap": true,
            "skip_rseqc": true,
            "skip_deseq2_qc": true,
            "extra_star_align_args": "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA",
            "extra_salmon_quant_args": "--noLengthCorrection"
        }
        ```
        
    - If including UMI deduplication:
        
        ```bash
        {
            "input": "<samplesheet.csv>",
            "outdir": "<output_dir>",
            "fasta": "<path_to_fasta>",
            "gtf": "<path_to_gtf>",
            "star_index": "<path_to_star_index>",
            "with_umi": true,
            "umitools_extract_method": "regex",
            "umitools_bc_pattern": "^(?P<umi_1>.{6})(?P<discard_1>.{4}).*",
            "umitools_dedup_stats": true,
            "skip_pseudo_alignment": true,
            "skip_markduplicates": true,
            "skip_bigwig": true,
            "skip_stringtie": true,
            "skip_dupradar": true,
            "skip_qualimap": true,
            "skip_rseqc": true,
            "skip_deseq2_qc": true,
            "extra_star_align_args": "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA",
            "extra_salmon_quant_args": "--noLengthCorrection"
        }
        ```
        
3. Create a run nextflow `.sh`  file (must be executable):
    
    ```bash
    #!/usr/bin/env bash
    workflow=/local/cqls/software/nextflow/assets/nf-core-rnaseq_3.22.2/3_22_2
    
    NXF_TEMP=/scratch nextflow run \
            $workflow \
            -work-dir <work_dir> \
            -profile singularity \
            -params-file <parameters.json> \
    ```
    
    - `<work_dir>` is where intermediate files will be stored
4. Ensure that `nextflow.config` is in the current working directory, then submit the run as a batch job by running `sbatch ./<run_nextflow>.sh`

## Results

### Monitoring Run

- To ensure that the run has submitted and has started without any errors, run `hqstat`
    - Initially, the only visible job running will be the shell script on `all.q`
    - After a few minutes, if the run has been configured correctly, more jobs will pop up on `biomed@ruby`, corresponding to individual processes within the Nextflow pipeline
- To monitor the run progress, view the contents of  `slurm-<runid>.out` and `<work_dir>`
    - If there are any errors, the `.out` log file will contain those
    - Run `tail -200 slurm-<runid>.out` to see the most recent progress
- Once the run has completed successfully, it will not appear in `hqstat` and the output log file will have the following lines at the end (exact numbers will be run-specific):
    
    ```bash
    Completed at: 15-Dec-2025 08:34:15
    Duration    : 5h 2m
    CPU hours   : 113.2
    Succeeded   : 970
    ```
    

### Pipeline Output

- All output files will be contained within `<output_dir>`
- Directories within `<output_dir>`:
    
    ```bash
    fastqc/
    fq_lint/
    multiqc/
    pipeline_info/
    star_salmon/
    trimgalore/
    
    # if included UMI deduplication
    umitools/
    ```
    
- The final gene count table will be contained within `star_salmon/`: `salmon.merged.gene_counts.tsv`
- Explore other directories for quality, trimming, and mapping statistics