# SHC 2021
# These commands were used to map expression to the transcriptomes described in 
#   Church and Extavour 2021
# The steps follow the agalma pipeline
# Additionally genetrees are annotated with Phyldog
# These commands were performed on the high performance cluster at Harvard University


# Map expression reads

#!/bin/bash
#SBATCH -n 24                   
#SBATCH -N 1                    
#SBATCH -t 7-00:00              
#SBATCH -p shared       
#SBATCH --mem=102G              
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=church@g.harvard.edu

cd /n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021
source activate agalma_Sep2018
export AGALMA_DB=$PWD/agalma/data/agalma.sqlite

cd agalma/scratch

agalma -t 24 -m 100G expression --id 008Ab 008D
agalma -t 24 -m 100G expression --id 012Ab 008D
agalma -t 24 -m 100G expression --id 024Bb 008D
agalma -t 24 -m 100G expression --id 008An 008D
agalma -t 24 -m 100G expression --id 012An 008D
agalma -t 24 -m 100G expression --id 024Bn 008D
agalma -t 24 -m 100G expression --id 008Ao 008D
agalma -t 24 -m 100G expression --id 012Ao 008D
agalma -t 24 -m 100G expression --id 024Bo 008D
agalma -t 24 -m 100G expression --id 016Bb 025A
agalma -t 24 -m 100G expression --id 016Cb 025A
agalma -t 24 -m 100G expression --id 028Ab 025A
agalma -t 24 -m 100G expression --id 016Bn 025A
agalma -t 24 -m 100G expression --id 016Cn 025A
agalma -t 24 -m 100G expression --id 028An 025A
agalma -t 24 -m 100G expression --id 016Bo 025A
agalma -t 24 -m 100G expression --id 016Co 025A
agalma -t 24 -m 100G expression --id 028Ao 025A
agalma -t 24 -m 100G expression --id 106Bb 106A
agalma -t 24 -m 100G expression --id 7_1_4 106A
agalma -t 24 -m 100G expression --id 7_2_1 106A
agalma -t 24 -m 100G expression --id 106Bo 106A
agalma -t 24 -m 100G expression --id 7_2_2 106A
agalma -t 24 -m 100G expression --id 106Bn 106A
agalma -t 24 -m 100G expression --id 7_1_1 106A
agalma -t 24 -m 100G expression --id 7_2_4 106A
agalma -t 24 -m 100G expression --id 8_0_15 055A
agalma -t 24 -m 100G expression --id 8_0_4 055A
agalma -t 24 -m 100G expression --id 8_0_9 055A
agalma -t 24 -m 100G expression --id 8_0_12 055A
agalma -t 24 -m 100G expression --id 8_0_2 055A
agalma -t 24 -m 100G expression --id 8_0_2_resequence1 055A
agalma -t 24 -m 100G expression --id 8_0_7 055A
agalma -t 24 -m 100G expression --id 8_0_16 055A
agalma -t 24 -m 100G expression --id 8_0_6 055A
agalma -t 24 -m 100G expression --id 040Bb 040C
agalma -t 24 -m 100G expression --id 040Bn 040C
agalma -t 24 -m 100G expression --id 040Bo 040C
agalma -t 24 -m 100G expression --id 040Db 040C
agalma -t 24 -m 100G expression --id 040Dn 040C
agalma -t 24 -m 100G expression --id 040Do 040C
agalma -t 24 -m 100G expression --id 18_0_14 040C
agalma -t 24 -m 100G expression --id 18_0_14_resequence1 040C
agalma -t 24 -m 100G expression --id 18_0_17 040C
agalma -t 24 -m 100G expression --id 18_0_17_resequence1 040C
agalma -t 24 -m 100G expression --id 16_1_4 16_1
agalma -t 24 -m 100G expression --id 16_2_4 16_1
agalma -t 24 -m 100G expression --id 16_1_2 16_1
agalma -t 24 -m 100G expression --id 16_2_2 16_1
agalma -t 24 -m 100G expression --id 16_1_1 16_1
agalma -t 24 -m 100G expression --id 16_2_1 16_1
agalma -t 24 -m 100G expression --id 002Cb 002D
agalma -t 24 -m 100G expression --id 032Ab 002D
agalma -t 24 -m 100G expression --id 032Bb 002D
agalma -t 24 -m 100G expression --id 002Cn 002D
agalma -t 24 -m 100G expression --id 032An 002D
agalma -t 24 -m 100G expression --id 032Bn 002D
agalma -t 24 -m 100G expression --id 002Co 002D
agalma -t 24 -m 100G expression --id 032Ao 002D
agalma -t 24 -m 100G expression --id 032Bo 002D
agalma -t 24 -m 100G expression --id CFAb CFB
agalma -t 24 -m 100G expression --id CFBb CFB
agalma -t 24 -m 100G expression --id CFCb CFB
agalma -t 24 -m 100G expression --id CFAn CFB
agalma -t 24 -m 100G expression --id CFBn CFB
agalma -t 24 -m 100G expression --id CFCn CFB
agalma -t 24 -m 100G expression --id CFAo CFB
agalma -t 24 -m 100G expression --id CFBo CFB
agalma -t 24 -m 100G expression --id CFCo CFB
agalma -t 24 -m 100G expression --id 020Bb 020A
agalma -t 24 -m 100G expression --id 020Cb 020A
agalma -t 24 -m 100G expression --id 020Db 020A
agalma -t 24 -m 100G expression --id 020Bn 020A
agalma -t 24 -m 100G expression --id 020Cn 020A
agalma -t 24 -m 100G expression --id 020Dn 020A
agalma -t 24 -m 100G expression --id 020Bo 020A
agalma -t 24 -m 100G expression --id 020Co 020A
agalma -t 24 -m 100G expression --id 020Do 020A
agalma -t 24 -m 100G expression --id 088Ab 088B
agalma -t 24 -m 100G expression --id 088Bb 088B
agalma -t 24 -m 100G expression --id 088Cb 088B
agalma -t 24 -m 100G expression --id 088An 088B
agalma -t 24 -m 100G expression --id 088Bn 088B
agalma -t 24 -m 100G expression --id 088Cn 088B
agalma -t 24 -m 100G expression --id 088Ao 088B
agalma -t 24 -m 100G expression --id 088Bo 088B
agalma -t 24 -m 100G expression --id 088Co 088B
agalma -t 24 -m 100G expression --id 021Ab 029A
agalma -t 24 -m 100G expression --id 029Db 029A
agalma -t 24 -m 100G expression --id 021An 029A
agalma -t 24 -m 100G expression --id 029Dn 029A
agalma -t 24 -m 100G expression --id 021Ao 029A
agalma -t 24 -m 100G expression --id 029Do 029A
agalma -t 24 -m 100G expression --id 043Cb 043D
agalma -t 24 -m 100G expression --id 056Ab 043D
agalma -t 24 -m 100G expression --id 043Cn 043D
agalma -t 24 -m 100G expression --id 056An 043D
agalma -t 24 -m 100G expression --id 043Co 043D
agalma -t 24 -m 100G expression --id 056Ao 043D
agalma -t 24 -m 100G expression --id 40_2_1 043D
agalma -t 24 -m 100G expression --id 40_2_4 043D
agalma -t 24 -m 100G expression --id 40_2_4_resequence1 043D
agalma -t 24 -m 100G expression --id 40_2_2 043D
agalma -t 24 -m 100G expression --id 40_2_2_resequence1 043D
agalma -t 24 -m 100G expression --id 40_2_2_resequence2 043D

# Annotate genetree with Phyldog

#!/bin/bash
#SBATCH -n 4                   
#SBATCH -N 1                    
#SBATCH -t 7-00:00              
#SBATCH -p serial_requeue       
#SBATCH --mem=48G              
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=church@g.harvard.edu

cd /n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021
source activate agalma_Sep2018
export AGALMA_DB=$PWD/agalma/data/agalma.sqlite

WORKDIR=/n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021/phyldog_results
CODEDIR=/n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021/phyldog_code

####### HERE IS WHERE YOU NEED TO FIND THE CORRECT MULTALIGN NUMBER #######
MULTIDIR=/n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021/agalma/scratch/multalign-55
###########################################################################

mkdir -p $WORKDIR
mkdir -p $WORKDIR/links
mkdir -p $WORKDIR/OptionFiles
mkdir -p $WORKDIR/ResultFiles

mkdir -p $CODEDIR

####### HERE IS WHERE YOU NEED TO FIND THE CORRECT TREE NUMBER #######
cp /n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021/agalma/scratch/iqtree_59/June2021_iqtree59_concat.treefile $CODEDIR/species_tree.tre
###########################################################################

cp /n/holyscratch01/extavour_lab/shchurch/PHYLDOG/parse_links.py $CODEDIR/
cp /n/holyscratch01/extavour_lab/shchurch/PHYLDOG/scripts/prepareData.py $CODEDIR/


### using the multalign from BEFORE treeprune
#   *    50  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 homologize    holy7c18609.rc.fas.harvard.edu shchurch 2021-06-15T01:34:48.030796    
#   *    51  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 multalign     holy7c18609.rc.fas.harvard.edu shchurch 2021-06-15T06:19:55.671542    
#   *    52  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 genetree      holy7c18609.rc.fas.harvard.edu shchurch 2021-06-16T18:42:10.897741    
#   *    53  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 treeinform    holy7c18609.rc.fas.harvard.edu shchurch 2021-06-16T19:06:27.315845    
#   *    54  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 homologize    holy7c18609.rc.fas.harvard.edu shchurch 2021-06-16T19:26:56.029560    
#   *    55  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 multalign     holy7c18609.rc.fas.harvard.edu shchurch 2021-06-17T00:07:32.028903    
#   *    56  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 genetree      holy7c18609.rc.fas.harvard.edu shchurch 2021-06-18T12:23:27.777957    
#   *    57  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 treeprune     holy7c18609.rc.fas.harvard.edu shchurch 2021-06-18T12:46:43.222883    
#   *    58  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 multalign     holy7c18608.rc.fas.harvard.edu shchurch 2021-06-19T00:37:00.228364    
#   *    59  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 supermatrix   holy7c18608.rc.fas.harvard.edu shchurch 2021-06-19T13:53:57.410041    
#   *    60  assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 supermatrix   holy7c26601.rc.fas.harvard.edu shchurch 2021-06-19T20:13:29.171250

cd $MULTIDIR
mkdir -p fa-gb # move fa-gb files for phyldog into here
cd alignments
for f in *.fa-gb
do
 NAME=${f%.fa-gb}
 OUT="$WORKDIR/links/$f"
 cat $f | python $CODEDIR/parse_links.py >$OUT
 cp $f ../fa-gb
done


cd $WORKDIR/
cut -d: -f1 $WORKDIR/links/* | sort -u >SpeciesNames.txt

# inputFile.txt for prepareData.py
echo $MULTIDIR/fa-gb > inputFile.txt
echo DNA >> inputFile.txt
echo fasta >> inputFile.txt
echo no >> inputFile.txt
echo $WORKDIR/links >> inputFile.txt
echo $WORKDIR/OptionFiles >> inputFile.txt
echo $WORKDIR/ResultFiles >> inputFile.txt
echo yes >> inputFile.txt # Give starting species tree?
echo $CODEDIR/species_tree.tre >> inputFile.txt #absolute path to starting species tree
echo no >> inputFile.txt # Optimize species tree?
echo yes >> inputFile.txt # Optimize DL parameters?
echo branchwise >> inputFile.txt # How to optimize DL parameters?
echo no >> inputFile.txt # Want to assume all species have same number of genes?
echo no >> inputFile.txt # Optimize gene trees?
echo yes >> inputFile.txt
echo 72 >> inputFile.txt # Time limit in hours

python2 $CODEDIR/prepareData.py < inputFile.txt

cd OptionFiles
for f in *.opt
do
 sed -i ':a;N;$!ba;s/input.sequence.sites_to_use=all/output.events.file=$(RESULT)$(DATA)_Events.txt \ninput.sequence.sites_to_use=all/' $f
done

cd $WORKDIR/

module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 boost/1.63.0-fasrc03
singularity exec ~/phyldog.simg mpirun -np 8 phyldog param=OptionFiles/GeneralOptions.txt


# Export expression

cd /n/holyscratch01/extavour_lab/shchurch/assembly_phylogeny_expression_Hawaiian_Drosophila_June2021
source activate agalma_Sep2018
export AGALMA_DB=$PWD/agalma/data/agalma.sqlite

agalma -t 1 -m 10G export_expression --id assembly_phylogeny_expression_Hawaiian_Drosophila_June2021 --speciestree phyldog_results/ResultFiles/StartingTree.tree --speciestree_numbered phyldog_results/ResultFiles/OutputSpeciesTree_ConsensusNumbered.tree --genetrees phyldog_results/ResultFiles/*.ReconciledTree --sequences 13 14 15 16 17 18 20 21 22 23 24 25  >agalma/reports/June2021_expression_export.json
