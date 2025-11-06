# nuclei

This is a work-in-progress guide on replicating my pupa developmental snRNASeq analysis. I'm trying to make it as painless to run as is possible, but I'm sure I've forgotten several key points. Let me know if you run into any issues you can't fix on your own.

Also, I'll note that you shouldn't just blindly run stuff. Read the code beforehand and try to vaguely understand what's going on. Particularly, I've hardcoded a lot of file locations that may or may not apply to you. Also, a lot of my scripts and the tools they call have fun options that I probably failed to mention here.

## Setup

I'm assuming you are running this on the lab workstation. With a little effort, everything here should be replicable on other Windows and/or Linux systems (Macs would take some work) with the appropriate hardware, i.e a lot of memory (maybe ~32GB, likely more) and a recent Nvidia (or AMD, with some hacks) dedicated GPU.

1. Install the latest versions of R, RStudio, and Bioconductor.
2. Install Ubuntu-24.04 under WSL. Instructions are [here](https://learn.microsoft.com/en-us/windows/wsl/basic-commands). It's important that you use 24.04 specifically and its possible that another version is already installed, so verify you have the correct one.
4. Download this repository. If you are comfortable using git, you can make a fork and use it to track any changes you make, or you can just download it as a zip file. I recommend placing it in the C or Z drives since they are much faster and you're about to be manipulating some large files. Drives D and Y are better for archival purposes. Unless stated otherwise, all scripts I tell you to run should be ran from `working_dir`.
5. In Ubuntu, install a python package manager (unless it already comes with one, I forget). My scripts are tested with [this one](https://conda-forge.org/download/), but others might also work.
6. Run create_cellbender_env.sh to install CellBender into a new conda environment called `cellbender`. Note that this installs my own customized version of CellBender which includes a bug fix for compatibility with our GPU, so if it breaks that's my fault.
7. TODO: whatever I forgot

## Mapping & Quantification

I already did this part for you, but here are some notes if you are interested.

The raw 10x 3'v3.1 reads are stored in `D:\Bai_lab_RNASeq_data\2024_Pupa_snRNAseq` and `D:\Bai_lab_RNASeq_data\Ping-White Pupa-Fasq`.

Each white pupa library was sequenced on two lanes. For the other stages, each library was sequenced 5 times with differing lane splits, read lengths, and read depths. These were all concatenated to yield one fastq pair per library.

I've experimented with several mappers: CellRanger, STARsolo, and alevin-fry. For now, I have settled on STARsolo, mostly because its the fastest (it only takes a week to run) and it exposes some advanced functionality that the others don't. However, upon more testing, I may end up switching to alevin-fry because it has some interesting UMI deduplication strategies. 

I aligned to the FlyBase v6.62 reference. 10x provides filtered genome annotations that may increase the unique mapping rate, but only for human and mouse. I'm working on a modified drosophila reference, but its going very slowly and it looks like I'll have to manually edit ~1000 gene records.

My STARsolo scripts are in `mapping_and_quantification`. You'll notice that I'm using a whole lot of nonstandard options. Most of these are to mimic the output of CellRanger or are settings for 3'v3.1 reads, but let me explain some of the others:

Basically, these are optimizations for future analysis of transposable elements. They allow more multimappers to be output to the BAM file, which you'll note is actually not being output. This is to save space because this would take ~2.5TB of storage.
```
--winAnchorMultimapNmax 100
--outFilterMultimapNmax 100
--outMultimapperOrder Random
```

This one uses expectation maximization to distribute counts from multimappers to the relevant genes. Normally, these are discarded, but some of the pupa stages have A LOT of multimappers, so I have been trying to retain them. This may or may not actually be a good idea.
```
--soloMultiMappers EM
```

The STARsolo output is in `Z:\Caelen\snRNAseq_v2`. If you look in each sample's Solo.out directory, you'll see several subdirectories. `Gene`, `GeneFull`, and `GeneFull_Ex50pAS` each have the quantification results. `GeneFull_Ex50pAS` is the one we care about. Read the STAR documentation if you want to know why. There is also `SJ` which has splice junction information if that's something you're interested in. Finally, we have `Velocyto`. This one splits the counts into categories based on whether or not they appear to have originated from spliced vs unspliced mRNA. This is useful for velocity analysis and some cell calling techniques.

In 'GeneFull_Ex50pAS', we have filtered counts, which have already had rudimentary empty droplet calling applied, and 3 types of raw counts. `raw` doesn't include multimappers, `raw_multi` does, and `raw_multi_floor` also does but the counts are guaranteed to be integers. I made this one since many downstream programs break with the non-integer counts introduced by EM.

I have also tabulated the various statistics output from all the STARsolo runs in `Z:\Caelen\snRNAseq_v2\STAR_stats.xlsx`. If any of you come up with any theories that explain the graphs in that file, I'm all ears.

## Empty Droplet Calling & Background Removal

So STAR does its own empy droplet filtering but it's 2025, we can use machine learning to do it. That's what CellBender does, while also removing background RNA for improved clustering and differential expression power.

I've included 2 different methods to run it, in scripts `run_cellbender_auto.sh` and `run_cellbender_expect_10k_tuned.sh`. CellBender has several parameters that can be tuned, and the auto script does what you'd expect. But, since P1 through P4 have many low count nuclei that are almost indistinguishable from the background, these are thrown out. If you'd like to keep these, `run_cellbender_expect_10k_tuned.sh` forces CellBender to expect on the order of 10000 non-empty nuclei, which is something it vehemently disagrees with, so you also have to tune some other parameters. Note that my tuned parameters may not work with this newest sequencing run, and they might not work with mutants. If you go down this route, you will likely have to retune them by hand, and possibly individually per sample.

For convieniently running CellBender on many samples, I have also included `run_cellbender_batch.sh`.

Once its done running, you will have the output in the working directory. Each one will have a report which must be closely scrutinized before continuing. Check out the [documentation](https://cellbender.readthedocs.io/en/latest/index.html) to see what you are looking for and how to fix it, if needed.

## Cell Calling