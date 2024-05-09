This is for notes!
## Git
### github push changes

cd Jinhui_Botany_563  
Nano/open [notebook-log.md](notebook-log.md)  
git add .  
git commit -m "informative message" # “name” can be changed according to need  
Git push 

### token 
Setting-developer setting-generate new token

## Dataset
### Filter & background
* Na+/K+-ATPase gene family sequence from three clade of Eurytemora affinis species complex
* Na+/K+-ATPase gene family sequence from other Arthropod (Chelicerata/Myriapoda/crustacea/hexapoda)  
* Ref sequences filter based on:  
  MINIMUM PROTEIN SIZE=200  
  MINIMUM COVERAGE=50  
  MINIMUM SCAFOLD N50=10000  
  MAX SCAFOLD COUNT=9000
* Lee lab has found signatures of parallel selection acting on the paralogs of NKA gene family during independent invasion from saltwater to freshwater.  
  The goal is to construct the evolutionary history and molecular evolution of this crucial ion transporter in Arthropoda. 
### Data
* E. affinis: 5 paralogs in Red clade; 5 paralog in Gulf clade; 6 paralog in Purple clade (1 of red clade and 1 of gulf clade is short and this may due to bad assembly/pseudogene)
* Arthropoda: 50 species across Chelicerata/Crustacea/Hexapoda/myriapoda; 129 paralogs and orthologues within Arthropoda

## Paralogs identification & BlastP

## Orthologs identification & Orthofinder
### local blastp
* previous I have used local blast to identify the paralogs and othologs based on similarity (>50) and e-value (0)  
* This is time consuming because need build every database for every species.
### Orthofinder
* Having compared with local blastp, this method is super fast and accurate to identify all paralogs and othologs (Diamond)
* The command line is also easy  



## Multiple sequence alignments (MSA) & T-coffee
https://tcoffee.org/Projects/tcoffee/documentation/index.html#  
(This T-coffee document provide detailed info about function and setting)
* T-Coffee offers a unique approach to multiple sequence alignment (MSA) by efficiently combining local and global alignment information. This capability makes it particularly effective for handling challenging cases, such as splice variants that require long gaps
* It provides versatility in handling different types of sequence data from DNA seq to protein seq. Also could provide different setting to diff seq type.
* Choosing the right package: T-Coffee is probably the most versatile, but it comes at a price, its default aligner being currently slower than many alternative packages. We can select which one to use according to our own usage.
* ! M-coffee: need more time to try.
* ! default setting: need more time to reading.

```
t_coffee your_sequences.fasta -type=dna -output=clustalw
```
PRACTICE: This is much more accurate than MAFFT!


## Phylogenetic tree
## Distance and parsimony methods for constructing tree using R



|             Software             |                                                                   Description                                                                    |                                                                     Strengths                                                                     |                                 Weaknesses                                  |                                                                                                                                   Assumptions                                                                                                                                   |                       User choices                       | 
|:--------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------:| 
|     R package: ape/adegenet      |                                                 non-Character-based data; <br/> distanced-based                                                  | get the tree in a scalable manner regardless of sample size;<br/>don't have to search the space of tree; <br/> fast;<br/> character-based methods | don't have to search the space of tree;<br/> depend on which model chosen   | Molecular Clock Assumption (mutations accumulate at a constant rate time/lineage); <br/> Additive Distances: (distance between two taxa (or sequences) is equal to the sum of the branch lengths along the path connecting them in the tree); <br/> Correct Distance Estimation | NJ; ME (minimum evolution); WPGMA; UPGMA Fitch algorithm |
| R package: phangorn/ape/adegenet | Character-based data; Parsimony-based method; <br/> seeks the tree that minimizes the amount of evolutionary change required to explain the data |                                         model-based method; <br/> look for trees; <br/> less computations                                         |                      ignore the real evolution process                      |                                           independence among characters；<br/> Ockham’s razor (the simple one to explain when equal explanations)                                                                                                                                |                    Fitch algorithm                       |


r for distance-based method
```r 
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
library(ape)
library(phangorn)

# This is distance-based method
# Ape is better for DNA tree constrcted
alignment <- read.phyDat("/Users/jinhuiwang/Desktop/20240305_tree_with_myripoda_v5.fasta", format = "fasta", type = "AA")
distMatrix <- dist.ml(alignment)
tree <- nj(distMatrix)
plot(tree, type="unrooted")
```
r for parsimony-based method
```r

# This is  Parsimony-based method
library(phangorn)
sequences <- read.phyDat("/Users/jinhuiwang/Desktop/20240305_tree_with_myripoda_v5.fasta", format = "fasta", type="AA")
#a starting nj tree for the search on tree space and compute the parsimony score of this tree
dist_matrix <- dist.ml(sequences)
nj_tree <- nj(dist_matrix)
# searching and modifying the Maximum Parsimony Tree
optimized_tree <- pratchet(sequences, nj_tree)
plot(optimized_tree, type = "unrooted")
tiplabels()
# output treefile
write.tree(optimized_tree, file="/Users/jinhuiwang/Desktop/tree_file.nwk")
```

## Maximum Likelyhood phylogenies
* How to improve current tree (uphill):  
  NNI nearest neighbor interchange (local optimal)  
  SPR subtree pruning and regrafting (great)  
  TBR tree bisection and reconnection (high computational demand)  
  Stochastic (随机) algorithms (allow "downhill")


|       Software        |                                                                                                                                                                                                                          Description                                                                                                                                                                                                                          |                                                                                                                                      Strengths                                                                                                                                      |                                                       Weaknesses                                                       | Assumptions |                                                                                                                                                    User choices                                                                                                                                                    | Reference                                                                                                                                   | 
|:---------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------:|:-----------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------------------------------| 
| IQ-TREE/IQ-TREE2<br/> | implementation of hill-climbing NNI (the possibility to find new and higher local optima) and stochastic NNI (random perturbation of locally optimal trees) operations; <br/> The extension of tree-search: maintain a small population (candidate tree set) of locally optimal trees initially generated from a large number of maximum parsimony (MP) trees. Throughout the ML tree search, we continuously update the candidate tree set with better trees |                                                                                    high efficient likelihoods than (RAxML and PhyML); <br/> more accuracy ultrafast bootstrap approximation<br/>                                                                                    | longer CPU time (for DNA) but have fast options; <br/> may finish at local optima and thus should rerun multiple times |             | -m "MIX{m1,...,mK}"  Mixture model with K components ; <br/> -mlrate site-specif rates ; <br/> -root-seq FILE,SEQ_NAME  Specify the root sequence from an alignment/non time-reversible evolutionary models (rooted tree); <br/> fast tree search to decrease time; -fast  Fast search (tree) to resemble FastTree | https://academic.oup.com/mbe/article/37/5/1530/5721363?login=false; https://academic.oup.com/mbe/article/32/1/268/2925592?login=false<br/>  |
|       RAxML-NG        |                                                                                                                                                                                               Evolutionary model extensions and chosen: 22 GTR-derived models                                                                                                                                                                                                 | Provides scalability and flexibility/ generally faster/ scores higher when alignment are taxon-rich/ use local tree arrangement/Multi-state morphological DNA/ RNA secondary structural DNA/ default is Maximum Parsimony but cal also random statrting tree and user provided tree |                                                                                                                        |             |                                                                                                                   cannot automaticly select best model/ (? only General time-reversible model )                                                                                                                    |                                                                                                                                             |




### Raxml-ng 
```
raxml-ng --check --msa 20240307_tree_with_myriapoda_v6.aln.fasta --model LG4M+G4 --log debug
raxml-ng --all --msa 20240307_tree_with_myriapoda_v6.aln.fasta --model PROTGTR --prefix T1 --seed 2 --threads 4 --bs-metric fbp,tbe 
```
* (specifies the number of starting trees --tree pars{N},rand{M})  
* https://isu-molphyl.github.io/EEOB563/computer_labs/lab4/raxml-ng.html  
* https://github.com/amkozlov/raxml-ng/wiki/Input-data#starting-trees

### iq-tree
```
iqtree -s new-hyphy-translatorx-manual2.nex -m TEST -bb 100000 -pre new-hyphy-copepoda-translatorx-manual2 -T 4
```
* try different times, they may have different results.

## Bayesian phylogenies
RevBayes & MrBayes
we could start from Mrbayes and then use the parameters in the Mrbayes to customize the parameters in the RevBayes.



## Bayesian analysis
### PHYLOBAYES-MPI
https://github.com/bayesiancook/pbmpi/blob/master/pb_mpiManual1.9.pdf
```
mpirun -np <n> pb-mpi -d datafile.ali  -cat -gtr -dgam 4 output2
nohup mpirun -np <n1> pb-mpi -d datafile.ali -cat -gtr -dgam 4 output_chain1 > output_chain1.log 2>&1 & nohup mpirun -np <n2> pb-mpi -d datafile.ali -cat -gtr -dgam 4 output_chain2 > output_chain2.log 2>&1 &
```
PAY ATTENTION  
at least 2 mpi: echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope  
2>&1： record both standard error and standard error  
  
**Advantages**  
Modeling  pattern-heterogeneity  across  sites variation in the rate of evolution (rate-heterogeneity across 
sites): diﬀerent rates & diﬀerent preferences for nucleotides or amino-acids (pattern-heterogeneity 
across sites), especially for long-branch attraction  

**global exchange rates setting**
CAT-Poisson/CAT-GTR/CAT  

**An initial tree can be provided**
input format:
PHYLIP formal (pay attention to the name)  

**Checking convergence and mixing**  
couldn't rely on the absolute number of cycles  
  
* _tracecomp program (for checking convergence of the continuous parameters of the model)_
```
tracecomp -x 1000 <chain1> <chain2> #.trace
```
chain 1 and chain 2 means run the program twice independently to make sure they began at the random seeds  
**ESS (effective sample size)**: independent sample and every sampled tree can provide riched info  
**rel diff** < 0.1 and minimum effective size > 300: good run  
**rel diff** < 0.3 and minimum effective size > 50: acceptable run



* _bpcomp program (for assessing convergence in tree space)_
```
bpcomp -x 1000 10 <chain1> <chain2> # using a burn-in of 1000, and sub-sampling every 10 trees # .treelist
```
**burn-in**: throw the data that is not real distribution on posterior distribution  
**maxdiff_ < 0.1**: good run  
**maxdiff < 0.3**: acceptable  
**0.3 < maxdiff < 1**: the sample is not yet sufficiently large and the chain have not convergaed, but this is on the right track   
if maxdiff =1 even after 10,000 points, this indicateds that at least one of the runs is stuck in a local maximum.  
**after bpcomp, you will get a tree document which include the _current best tree_

* _others_  
-t    <treefile>  
forces  the  chain  to  start  from  the  speciﬁed  tree  
-T    <treefile>  
forces the chain to run under a ﬁxed topology (as speciﬁed in the given ﬁle). In other words, 
the chain only samples from the posterior distribution over all other parameters (branch 
lengths, alpha parameter, etc.), conditional on the speciﬁed topology. This should be a 
bifurcating tree (see Input Data Format section).


## Mrbayes (not compatible with mac M1)
**links**  
https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf  
https://biojuse.com/2023/07/06/%E8%B4%9D%E5%8F%B6%E6%96%AF%E5%BB%BA%E6%A0%91%E4%B9%8B%20Mrbayes%20%E7%AF%87/


**_Input data_**  
Nexus (uses fixed set of symbols for each type of data)

**_specifying the model (e.g. GTR+I+T)_**  
_IQTREE ModelFinder to select the best model for your data_

```
MrBayes > log start filename = mrbayes.log # preserve screen output to mrbayes.log
> Mrbayes execute <filename.nex> # in same folder or absolute path
> Mrbayes help lset # check help and after setting
> Mrbayes lset nst=6 rates=invgamma
```  

_**setting the Priors (work well for most analysis)**_  
_add the mrbayes block to the end of the nex file_
```
> Mrbayes prset 
```
Topologypr: the topology  
Brlenspr: the branch length  
Statefreqpr: the four stationary frequencies of nucleotides  
Revmatpr: the six diﬀerent nucleotide substitution rates  
Pinvarpr: the proportion of invariable sites  
Shapepr: the shape parameter of the gamma distribution of rate variation.  

_OR the default preset_

```
begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=10000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
 outgroup Anacystis_nidulans;
 mcmc;
 sumt;
end;
```

**_Checking the Model_**
```
> Mrbayes showmodel
```

**_setting up the analysis_**  
_the chains is relavent to the nchains (e.g. under 4 default chains, 3 cold chain and 1 hot chain) and nruns;_  
_[cold chain] (hot chain)_

```
> Mrbayes help mcmc
> Mrbayes mcmcp ngen=20000 nchains=10 Nruns=2 diagnfreq=1000 samplefreq=100 printfreq=100
```

```
0 -- [-8849.255] (-8762.724) (-8987.103) (-8764.951) * [-8802.317] ...
```
**_running the analysis_**
```
> Mrbayes mcmc
```

**_when to stop/convergence  diagnostic_**  
* _ASDSF examine the average standard deviation of split frequencies_  
< 0.01: very good indication of convergence  
0.01-0.05: may be adequate depending on the purpose of your analysis  
* OR examine the log likelihood values (or, more exactly, the log probability of the data given the parameter values) of the cold chain  
* OR the mixing and convergence of the MCMC
PSRF Potential scale reduction factor: examine the different chain the intra-chain variance and the inter-chain variance whether are on the same distribution  
1-1.2  
ESS > 200 is good
```
> Mrbayes sump .p
```

**_summerizing the tree sample .t_**  
```
MrBayes > sumt conformat=Simple contype=Halfcompat relburnin=yes burninfrac=0.25
```  
contype=halfcompat/Allcompat
sumt: allow to examine the error/standard deciation associated with each clade in the tree (less supported clade: posterior probabilities <0.95)

## The coalescent method-ASTRAL-3
https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md  
https://arxiv.org/abs/1904.03826
* install check
```
 java -jar astral.5.7.8.jar -i test_data/song_primates.424.gene.tre
```
```
java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
java -Djava.library.path=./lib/ -jar astralmp.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
```
! need more data to better infer the species tree



## view tree software  
iTOL/ggtree/figtree/Treeviewer
## change tip labels 
https://stackoverflow.com/questions/57547769/change-tip-labels-in-phylogeny-in-r  
```
library("ape")

orig_tiplabels <- c("Alice", "Bob", "Cindy")
orig_tree <- rtree(n = 3, tip.label = orig_tiplabels)
plot(orig_tree)

new_tiplabels <- c("Debbie", "Elrond", "Frank")
orig_tree$tip.label <- new_tiplabels
plot(orig_tree)
```

## Others
* estimate time/sepciation age: need to estimate the time after infer the tree
* Phylogeny tree can be affected by a lot of reasons. The data is not just sequences: must be accurate and good quality sequences. But some protein is diverse. **_SO, YOU COULD LEAVE IDENTITY > 70 SEQUENCES RO MAKE SURE GOOD ALIGHMENT AND COULD REFLECT THEIR TRUE RELATIONSHIP._**