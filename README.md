# MCTS-RNA
MCTS-RNA is a tool for RNA inverse folding problem based on Monte Carlo Tree Search method. MCTS-RNA can design nested RNA structures and pseudoknot structures with user designed constraints: wide range and precise GC-content constraint and GC-content devation constraint. 
# Requirements
1.RNAfold of [ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/index.html) needs to be installed for designing nested structures.

2.[pKiss](http://bibiserv2.cebitec.uni-bielefeld.de/pkiss) is required for designing pseudoknot structures. 

2.[Python](https://www.python.org/) required version is at least version 2.7.

3.Python library of [Numpy](http://www.numpy.org/) need to be installed.

# Installation
You can download the python script MCTS-RNA.py, run this script from the shell. 


# How to use MCTS-RNA?
Once you downloaded the python script of MCTS-RNA and having installed all the requirements, you can execute MCTS-RNA from the shell. The inputs include the dot-bracket representation of target RNA secondary structure ,the target GC-content of the RNA sequence and GC-content error. The following are the examples and explanations of the inputs parameters.

This is an example of the command in the shell for nested RNA structures.

python MCTS-RNA.py -s "...(((((..........)))))........((((((((......))))))))(((((.......))))).............(((((..(((((..((..((.(((((.(((((.......))))).)))))...))....))))))))))))" -GC 0.75 -d 0.01 -pk 0

This is an example of the command in the shell for pseudoknot RNA structures.

python MCTS-RNA.py -s "....(((((.[[[[.))))).........]]]]..." -GC 0.4 -d 0.02 -pk 1


-s : The target RNA secondary structure.

-GC: The target GC-content of the RNA sequence, choose vaule from the range [0,1]. 

-d : The GC-content deviation of the solution, which is in range [0,0.02]. MCTS-RNA can output the sequence with more accurate GC-content with smaller GC-content devation setting, the default value of the GC-content devation is 0.01.

-pk: Design nested structure by setting -pk 0 and design pseudoknot structure by setting -pk 1 (currently MCTS-RNA only uses Pkiss) to predict pseudoknot structures. 

