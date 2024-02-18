This is for notes!
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



## Git
### github push changes

Cd Jinhui_Botany_563 
Nano/open [notebook-log.md](notebook-log.md)  
git add .  
git commit -m "informative message" # “name” can be changed according to need  
Git push 

### token 
Setting-developer setting-generate new token