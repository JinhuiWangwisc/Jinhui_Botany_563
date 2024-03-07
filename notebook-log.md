This is for notes!
## Git
### github push changes

Cd Jinhui_Botany_563 
Nano/open [notebook-log.md](notebook-log.md)  
git add .  
git commit -m "informative message" # “name” can be changed according to need  
Git push 

### token 
Setting-developer setting-generate new token

## Dataset
* Na+/K+-ATPase gene family sequence from three clade of Eurytemora affinis species complex
* Na+/K+-ATPase gene family sequence from other Arthropod (Chelicerata/Myriapoda/crustacea/hexapoda)  
* Ref sequences filter based on:  
  MINIMUM PROTEIN SIZE=200  
  MINIMUM COVERAGE=50  
  MINIMUM SCAFOLD N50=10000  
  MAX SCAFOLD COUNT=9000
* Lee lab has found signatures of parallel selection acting on the paralogs of NKA gene family during independent invasion from saltwater to freshwater.  
  The goal is to construct the evolutionary history and molecular evolution of this crucial ion transporter in Arthropoda. 

## Multiple sequence alignments (MSA) & T-coffee
https://tcoffee.org/Projects/tcoffee/documentation/index.html#  
(This T-coffee document provide detailed info about function and setting)
* T-Coffee offers a unique approach to multiple sequence alignment (MSA) by efficiently combining local and global alignment information. This capability makes it particularly effective for handling challenging cases, such as splice variants that require long gaps
* It provides versatility in handling different types of sequence data from DNA seq to protein seq. Also could provide different setting to diff seq type.
* Choosing the right package: T-Coffee is probably the most versatile, but it comes at a price, its default aligner being currently slower than many alternative packages. We can select which one to use according to our own usage.
* ! M-coffee: need more time to try.
* ! default setting: need more time to reading.

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
