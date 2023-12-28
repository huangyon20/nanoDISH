# nanoDISH
## Deconvolute the RNA intermediates structure heterogeneity with Nanopore direct RNA sequencing data.
Nanopore direct RNA sequencing that combined with chemical probing method was used to predict the secondary strucure of long RNA. This pepilne start with the raw fastq reads and fast5 signals of both modified and unmodified samples. Reads selection and One-class SVM were used to predict the modified bases of each selected reads that belong to a specific intermadiate. We transfer the modify profile of every selected reads to a bitvector and cluster the reads into different groups. The reactivity scores of each group, which represent the alternative comformations , were calculated and normalized individually. 

![flow](docs\Figures\Flow.png)

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
In this step, the reads in fastq format were mapped to reference genome using grphmap2 or minimap2. Then align the events from raw fast5 to the reference. The erroneously called redundant events were collapsed into one with recalculating the values of Current, Current sdtv and Dwelltime. The output is an collapsed-event file. The modified and unmodified group should be deal with separatly. The dependent softwares are:

```
python dataprocessing.py  -r ref.fasta -q mod.fastq   -s mod_fast5   -d Mod   -n mod   --mkindex
python dataprocessing.py  -r ref.fasta -q unmod.fastq -s unmod_fast5 -d Unmod -n Unmod --mkindex
```

step 2: Pick the intermidiates. 
-------------------------------
In this step, the reads in the cwere selected according to the mapping positions on refenrence. -i: the collapsed event file from the output of step 1. -s: the start mapping position. Reads were excluded if the start mapping postion is less than s. This limitation factor can exclude the degraded reads. -e: the end mapping position; -b: the allowed bias of end mapping positions. e.g if e=100 and b=20, which means the reads whose end mapping position ranged in 100-120 were specified to one intermediate.  If users want to select all the intermediates that shorter than the specific intermediate at the same time, -m should be used. e.g if s=30,e=100,b=20 and -m, which means the reads whose end mapping position ranged in (30-40,40-60,60-80,80-100,100-120) were specified to 5 intermediate individually.

```
python picksize.py -i Mod/mod.collapsed.event  -g TPP -d Mod/mod_multi -s 30 -e 340 -b 20 -m
python picksize.py -i Unmod/Unmod.collapsed.event  -g TPP -d Unmod/Unmod_multi -s 30 -e 340 -b 20 -m
```

step3: Predict modified bases.
--------------------------------------------
In this step, The intermediate event files from both modified and unmodified group were choosen to predict the modified bases using SVM. The output files (*.static and *.Mod_profile in csv format )were written to the path of --Mod_Event. If users wants to change the usage of SVM, there are 4 parameters available. --columns, 4: current, 5: current stdv, 6: dwell time, the default is [4,5,6]; --kernel, default is (rbf); --gamma, default is (scale); --nu, default is (0.01).

```
python SVM.py -m Mod/mod_multi/TPP_340.event -u Unmod/Unmod_multi/TPP_340.event -r ref.fasta -s TPP
```

step4: Calculate the reactivity scores.
--------------------------------------------
Based on the mod_profile file output in the step3, we calculate the modify-ratio for each base. Then use the Winsorization algorithm to normalize the ratios to reactivity scores along the full length of selected intermediate. The reactivity scores can be calculated generally as an average(m=mean), or as individual values taking into account alternative conformations(m=heter). 
```
python score.py -i Mod/mod_multi/TPP_340.event.Mod_Profile -o Mod/mod_multi/TPP_340_mean -r TPP.fa -m mean
python score.py -i Mod/mod_multi/TPP_340.event.Mod_Profile -o Mod/mod_multi/TPP_340_heter -r TPP.fa -m heter
```
During the running of the program, the user needs to view the temporarily produced PCA file and specify the number and method of clustering.

