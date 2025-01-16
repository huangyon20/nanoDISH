# nanoDISH
## Deconvolute the RNA intermediates' structure heterogeneity from Nanopore direct RNA sequencing data.
Nanopore direct RNA sequencing that integrated with chemical probing method was used to determine the secondary strucure of long RNA. Chemical probing reagents, such as NAI-N3 and Aclm, preferentially modified the open RNA bases and alter the current signals, resulting in miscalled bases or correctly called bases but with current signal outliers compared with unmodified bases. We used machine learning models to obtain the modification status of each bases on reads from the Nanopore sequencing data and developed this pipeline to analyze RNA structural heterogeneity.
NanoDISH consists of four parts: process the signal data, select reads by their length, obtain the modification profile of each modified reads, classify reads and calculate the reactivaty scores. The modification profiles were converted into 0-1-2 digital vectors and classified using Gaussian Mixture or DBSCAN method. The reactivity scores of each group were calculated and normalized individually. 

![flow](docs/Figures/Flow2.png)

### The dependent programs to run the analysis:
- graphmap2(https://github.com/lbcb-sci/graphmap2) or minimap2(https://github.com/lh3/minimap2) to align the Nanopore reads to reference.
- samtools(https://www.htslib.org) to deal with the sam/bam files.
- nanopolish(https://github.com/jts/nanopolish)


### The dependent Python packages:
- os
- argparse
- pandas
- numpy
- OneClassSVM in sklearn.svm
- PCA in sklearn.decomposition 
- GaussianMixture in sklearn.mixture
- DBSCAN in sklearn.cluster
- mean in statistics


step 1: Dataprocessing.
--------------------------------------------
In this step, the fastq reads were mapped to reference genome using grphmap2 or minimap2 and the fast5 signals were aligned to the reference RNA using nanopolish. The erroneously called redundant signals in event file were collapsed into one with recalculating the mean value of Current_mean, Current_sdtv and the sum value of Dwell time. The output file is an event file with no redundant signals. The modified and unmodified group should be deal with separatly. 

```
python dataprocessing.py  -r ref.fasta -q mod.fastq   -s mod_fast5   -d Mod   -n mod 
python dataprocessing.py  -r ref.fasta -q unmod.fastq -s unmod_fast5 -d Unmod -n Unmod

```

step 2: intermidiates selection. 
-------------------------------
If the raw sequencing reads are mixed with various length of RNA intermediates, this step could separate the reads into individual groups which corresponding to intermediates. The reads were grouped according to their mapping positions on refenrence. -i: the event file from step 1. -s: the start mapping site. Reads whose 5 terminal digested beyond the start site were excluded. -e: the end mapping site; -b: the allowed bias of end mapping positions. e.g., if e=100 and b=20, which means the reads whose end mapping sites ranged in 100-120 belong to one intermediate.  If users want to select all the intermediates that shorter than the specific intermediate at the same time, -m should be used. e.g., if s=30,e=100,b=20 and -m, which means the reads whose end mapping position ranged in (30-40,40-60,60-80,80-100,100-120) were specified to 5 intermediates individually.

```
python pickIntermediate.py -i Mod/mod.collapsed.event  -g TPP -d Mod/mod_multi -s 30 -e 340 -b 20 -m
python pickIntermediate.py -i Unmod/Unmod.collapsed.event  -g TPP -d Unmod/Unmod_multi -s 30 -e 340 -b 20 -m
```

step3: Making modification profiles.
--------------------------------------------
In this step, The event files from both modified and unmodified group were used to obtain the modification status of each bases. The output files (*.static and *.Mod_profile in csv format )were written to the path of --Mod_Event. The signal features used in the SVM model could be assigned though --columns (the default value is [4,5,6], where 4 represent for Current_mean, 5 for Current_stdv and 6 for Dwell time). The default --kernel is (rbf), the default --gamma is (0.01) and the default --nu is (0.01).

```
python predictModbase.py -m Mod/mod_multi/TPP_340.event -u Unmod/Unmod_multi/TPP_340.event -l 340 -b 20 
```

step4: Calculate the reactivity scores.
--------------------------------------------
We calculate the modify-ratio for each base based on the modification profile form the results of step3, then normalize the ratios to reactivity scores using Winsorization algorithm. The reactivity scores can be calculated generally as an average(m=mean), or as individual values (m=heter). 
```
python scoreIntermediate.py -i Mod/mod_multi/TPP_340.event.Mod_Profile -o Mod/mod_multi/TPP_340_mean -l 340 -b 20 -m mean
python scoreIntermediate.py -i Mod/mod_multi/TPP_340.event.Mod_Profile -o Mod/mod_multi/TPP_340_heter -l 340 -b 20  -m heter
```
During the running of the program, the user needs to view the temporarily produced UMAP file and specify the cluster number and clustering method.

