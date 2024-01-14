# Pipeline for bioinformatic analysis

Title - *A matter of salt: global assessment of the effect of salt ionic composition as a driver of bacterial diversity*\
Author - Attila Szabó\
Created On - 14 January 2024

## S3 Supplementary scripts of the bioinformatic analysis

Since sequence data were obtained by different sequencing platforms and chemistries, the versatile mothur tool (Schloss et al., 2009) was applied to process the raw sequences using a custom-created pipeline to obtain a unified high-quality OTU table using 99% similarity cutoff (see Supplementary Methods). 
These mothur scripts outline the pipeline for processing V4 and V3-V4 16S rRNA gene amplicon data using mothur *v1.44.3*. However, the preprocessing of 454 pyrosequencing reads was carried out with mothur *v1.35*.
This script also can be found on the GitHub repository: https://github.com/attiszabo/MatterOfSalt/

**Sequence set SRA accesions used in this study are listed in Supplementary Table 1.**


### Note
There are a few differences from mothur version 1.47 on the handling and creations of the name, count and group files

## Preprocessing

### Preprocess 454 pyrosequencing reads

#### Convert SFF to fasta and quality files

```{bash}
sffinfo(sff=run_name.sff)
```

#### Trim flowgrams

```{bash}
trim.flows(flow=run_name.flow, oligos=run_name_barcodes_341F_785R.oligos, minflows=350, maxflows=700, pdiffs=1)
  #515F_806R.oligos used for the V4 dataset, 341F_805R.oligos for the V3-V4 dataset
```

#### Denoise flowgrams

```{bash}
shhh.flows(file=run_name.flow.files)
```

#### Trim barcode, primers and quality filter
```{bash}
trim.seqs(fasta=current, name=current, oligos=run_name_barcodes_341F_785R.oligos, pdiffs=1, bdiffs=0, maxhomop=8, flip=T)
  #515F_806R.oligos used for the V4 dataset, 341F_805R.oligos for the V3-V4 dataset
```

*oligos file*(https://mothur.org/wiki/oligos_file/):\
barcode   *SAMPLE1_BARCODE*   SAMPLE1\
barcode   *SAMPLE2_BARCODE*   SAMPLE2\
barcode   *SAMPLE3_BARCODE*   SAMPLE3\
forward   CCTACGGGNGGCWGCAG   341F\
reverse   GACTACHVGGGTATCTAATCC  785R\


#### Summarize sequences
```{bash}
summary.seqs(fasta=current, name=current)
```

**Using mothur v1.44.3:**

#### Find and list unique sequences
```{bash}
unique.seqs(fasta=current)
```

#### Create a count_table
```{bash}
count.seqs(name=current, group=current)
```

#### Summarize sequence and group (sample) information
```{bash}
summary.seqs(fasta=current, count=current)
count.groups(count=current)
```

###  Preprocess pair-end Illumina sequences

#### List samples and sequence files
```{bash}
make.file(inputdir=., type=gz, prefix=project_name)
```

#### Combine paired reads of forward and reverse direction
```{bash}
make.contigs(file=project_name.files, deltaq=10)
```

#### Summarize sequences
```{bash}
summary.seqs(fasta=project_name.trim.contigs.fasta)
```

#### Count sequences per samples
```{bash}
count.groups(group=project_name.contigs.groups)
```

#### Quality filter reads
```{bash}
screen.seqs(fasta=current, group=current, maxambig=0, minlength=200, maxhomop=8)
	#'minlength' adjusted 200 for the V4 dataset, 400 for the V3-V4 dataset
```

#### Summarize filtered reads
```{bash}
summary.seqs(fasta=current)
```

### Preprocess non-paired-end Illumina sequences

#### Create a fasta and quality file from fastq
```{bash}
fastq.info(fastq=SRRxxxxxxx.fastq)
rename.file(fasta=current, qfile=current, prefix=sample1)
```

#### Create a group file and count reads per sample
```{bash}
make.group(fasta=sample1.fasta-sample2.fasta-sample3.fasta, groups=sample1-sample2-sample3)
count.groups(group=current)
```

#### Merge fasta files within the same project
```{bash}
merge.files(input=sample1.fasta-sample2.fasta-sample3.fasta, output=non_paired_project_name.fasta)
summary.seqs(fasta=current)
```

#### Quality filter reads
```{bash}
screen.seqs(fasta=current, group=current, maxhomop=8, maxambig=0, minlength=200, maxlength=300)
	#'minlength' adjusted to 200 for the V4 dataset, 400 for the V3-V4 dataset
	#'maxlength' adjusted to 300 for the V4 dataset, 500 for the V3-V4 dataset
summary.seqs()
```

#### Primer removal from Illumina datasets
```{bash}
trim.seqs(fasta=current, oligos=515F_806R.oligos, pdiffs=2, checkorient=T)
	#515F_806R.oligos used for the V4 dataset, 341F_805R.oligos for the V3-V4 dataset
summary.seqs()
# Optionally, if trim.seqs eliminated some of your sequences you need to list and remove them from your group file
list.seqs(fasta=project_name.trim.contigs.good.scrap.fasta)
remove.seqs(group=project_name.contigs.good.groups, accnos=project_name.trim.contigs.good.scrap)
# For cases where unidentified adapter residuals or primer sequences remained, mothur's trim.seqs was used in combination with the MEGA version 6 or 11 (Tamura, Stecher, Peterson, Filipski, and Kumar 2013) on previously SILVA 138 SSU aligned fasta files 
```
#### Find and list unique sequences
```{bash}
unique.seqs(fasta=current)
```

#### Create a count_table
```{bash}
count.seqs(name=current, group=current)
```
#### Summarize sequence and group (sample) information
```{bash}
summary.seqs(fasta=current, count=current)
count.groups(count=current)
```

# Merging datasets

#### For the V3-V4 dataset

```{bash}
merge.files(input=V3V4.project_name1.fasta-V3V4.project_name2.fasta-V3V4.project_name3.fasta…V3V4.project_namen.fasta, output=V3V4.fasta)
merge.count(count=V3V4.project_name1.count_table-V3V4.project_name2.count_table-V3V4.project_name3.count_table…V3V4.project_namen.count_table, output=V3V4.count_table)
summary.seqs(fasta=current, count=current)
count.groups(count=current)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)
```

#### For the V4 dataset

Ahead of merging, the MEGA software was used to trim the 341-533 region from the pre-aligned V3-V4 dataset, ensuring compatibility with the V4 dataset.

```{bash}
merge.files(input=V3V4_533F_trim.fasta-V4.project_name1.fasta-V4.project_name2.fasta…V4.project_namen.fasta, output=V4.fasta)
merge.count(count=V3V4.count_table-V4.project_name1.count_table-V4.project_name2.count_table…V4.project_namen.count_table, output=V4.count_table)
summary.seqs(fasta=current, count=current)
count.groups(count=current)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)
```

# Processing quality filtered and trimmed, merged datasets (V4/V3-V4)

#### Sequence Alignment

```{bash}
align.seqs(fasta=current, reference=/path/Silva_nr138/silva.nr_v138.align)
summary.seqs(fasta=current, count=current)
```

#### Removing non-overlapping reads
```{bash}
screen.seqs(fasta=current, count=current, start=13862, end=23440)
	#'start' adjusted to 13862 for the V4 dataset, 6428 for the V3-V4 dataset
summary.seqs(fasta=current, count=current)
count.groups(count=current)
```

#### Removing overhanging positions
```{bash}
filter.seqs(fasta=current, trump=., vertical=T)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)
```

#### De-noising reads
```{bash}
pre.cluster(fasta=current, count=current, diffs=2)
  #'diffs' adjusted to 2 for the V4 dataset, 4 for the V3-V4 dataset
summary.seqs(fasta=current, count=current)
count.grouos(count=current)
```

#### Chimera check and  removal
```{bash}
chimera.vsearch(fasta=current, count=current, dereplicate=T)
remove.seqs(fasta=current, accnos=current)
summary.seqs(fasta=current, count=current)
count.groups(count=current)
```

#### Singleton removal
```{bash}
split.abund(fasta=current, count=current, cutoff=1)
summary.seqs(fasta=\*.abund.fasta, count=\*.abund.count_table)
count.groups(count=current)
```

#### Taxonomic assignment and the removal of non-primer specific targets
```{bash}
classify.seqs(fasta=current, count=current, reference=silva.nr_v138.align, taxonomy=silva.nr_v138.tax, cutoff=80, method=wang)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Archaea-Chloroplast-Mitochondria-Eukaryota-unknown)
summary.seqs(fasta=current, count=current)
count.groups(count=current)
#For a more robust comparison, removing samples with less than 2500 reads
remove.groups(fasta=current, count=current, groups=sample_x-sample_y-sample_z)
count.groups(count=current)
```

# OTU-based analyses

#### Calculate distances and clustering sequences to OTUs
```{bash}
dist.seqs(fasta=current, cutoff=0.1)
cluster(column=current, count=current, cutoff=0.01)
```

#### Create OTU99 table and classify OTUs
```{bash}
make.shared(list=current, count=current, label=0.01)
classify.otu(list=current, label=0.01, taxonomy=current, count=current, basis=sequence)
```

#### Get OTU representative sequences (based on abundance)
```{bash}
get.oturep(list=current, fasta=current, count=current, method=abundance)
```

#### Calculate alpha diversity metrics
```{bash}
summary.single(shared=current, subsample=T, calc=nseqs-coverage-sobs-ace-chao-simpson-shannon-invsimpson-shannoneven-simpsoneven)
```

#### Rarefying OTU table for statistical analyses
```{bash}
sub.sample(shared=current)
```

# More-detailed OTU taxonomy using a freshwater database with SILVA 

Additional Taxonomic assignments of the OTU99 representatives for both datasets were made using TaxAss (Rohwer et al., 2018 Msphere) with the FreshTrain 2020Jun15 and Silva SSU v138 databases

- Using a text editor header of the \*.0.01.rep.fasta were modified with regular expressions to contain only the OTU IDs
- file was renamed to 'otus.fasta' as an input for TaxAss

```{bash}
cd TaxAss-master/tax-scripts/scripts
./RunSteps_quickie.sh otus FreshTrain15Jun2020silva138 silva_nr_v138_taxass 98 80 80 8
```
