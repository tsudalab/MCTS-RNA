# MCTS-RNA
MCTS-RNA is a tool for RNA inverse folding problem based on Monte Carlo Tree Search method. MCTS-RNA can satisfy with user designed constraints: wide range and precise GC-content constraint and GC-content error constraint. 
# Requirements
1.For the usage of MCTS-RNA, the program RNAfold of the ViennaRNA Package version 2.1.9 with python interface are required.
They need to be installed in the PATH variable of your computer.

2.Python required version is at least version 2.7

3.Python library of Numpy need to be installed.
#Installation
You can download the python script MCTS-RNA-Final.py and save it in the executable paths of your computer ,run this script from the shell. 
#How to use MCTS-RNA?
Once you downloaded the python script of MCTS-RNA and having installed all the requirements, you can execute MCTS-RNA from the shell. The inputs include the dot-bracket representation of target RNA secondary structure and the target GC-content of the RNA sequence, target GC-content and GC-content error.The following are the examples and explanations of the inputs parameters.

This is an example of the command in the shell:
python MCTS-RNA-Final.py -f -s "...(((((..........)))))........((((((((......))))))))(((((.......))))).............(((((..(((((..((..((.(((((.(((((.......))))).)))))...))....))))))))))))" -GC 0.75 -d 0.01

-s : The target RNA secondary structure.

-GC: The target GC-content of the RNA sequence, choose vaule from the range [0,1]. The default value of the target GC-content is random, which means the GC-content of the output solution is any value in range [0,1]. 

-d : The GC-content error of the solution, which is in range [0,0.02]. MCTS-RNA can output the sequence with more accurate GC-content with smaller GC-content error setting, the default value of the GC-content error is 0.02.
#Example of the output
1.Solution with the GC-content constraint and GC-content error constraint

python MCTS-RNA-Final.py -f -s "...(((((..........)))))........((((((((......))))))))(((((.......))))).............(((((..(((((..((..((.(((((.(((((.......))))).)))))...))....))))))))))))" -GC 0.75 -d 0.01

search length:112

Solution:GACGCGCCUUUAAGUUUUGGCGCCGCCUCUCGCGCCCCCCGUGCUGGGGGCGCAGGGCCGCAUACGCCCUCAACGCUAGAAUACGGCCUCGGGCCCAGCAACCCCCCGCCGGGGGCCCUGAACCCCCUGCGGGUCAGGAACGGCGGCCCGGCCG

running time:4.51990914345

GC content:0.746753246753

GC distance:0.00324675324675

structure distance:1.0


2.Solution without the GC-content constriant

python MCTS-RNA-Final.py -f -s "...(((((..........)))))........((((((((......))))))))(((((.......))))).............(((((..(((((..((..((.(((((.(((((.......))))).)))))...))....))))))))))))" 

search length:112

solution:UUCCGGGCCUCUUCGUAUGCCCGAAAUAAUUAUCCGUCCUCUAGUGGACGGAUCCGGGGACUUAACCCGGUCUUUCGCCGCUCGGUACUGGUCGGAUGGAAGCUUAUCUACCUUGGUCUACGCAAGGCAGAUACACGCGUAGCCCCGACGUACC

running time:0.27837896347

GC-content:0.564935064935

structure distance:1.0
