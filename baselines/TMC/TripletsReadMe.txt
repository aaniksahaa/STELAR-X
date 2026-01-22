Installation and running the program "treeFromTriplets" 

This ReadMe file describes the installation and the usage of the treeFromTriplets software. This software generates a phylogenetic tree using Weighted Rooted Triplet MaxCut Algorithm. The theoretical background for this algorithm is given in the article "Triplet MaxCut: A Fast Algorithm for Amalgamating Triplet Trees", authors: Zeev Frenkel and Sagi Snir, 2015.

This version is compiled under x86_64 GNU/Linux
Ubuntu 10.04.4 LTS

Input/Output formats:

0. Input tree(s) - Optional. Single file containing single rooted tree (or trees) in Newick format (see http://en.wikipedia.org/wiki/Newick_format). Each of leaves should be named and only leaves should be named. Distances are not play role (can be omitted). Each of input tree should be followed by symbol ";". There are no formal limitations on the maximal length of leaf name, number of trees and number of leaf names (depend on the system). Please note that in Newick format is not allowed usage of symbols "(", ")", ":", ";" and "," in leaf names. It is also not allowed to use symbol "," in the distances, i.e., non-integer distances should be written using decimal dot ".", not comma ",". Rooting of input tree(s) is considered as a default one for the Newick format, i.e., node with name "root", "out group" etc. is considered as regular leaf. 

1. Input triplets - Optional. Single file containing the set of rooted triplets in format: <index of leaf1>, <index of leaf2> | <index of leaf3>[:<weight>]. Indexes of leaves should be non-negative integer numbers corresponding to taxa. Weight is an optional parameter, positive number (non-integer weights should be written using decimal dot ".", not comma ","). Please don't write weights for triplets in the case of algorithm usage without weights. Please write weight for each rooted triplet in the case of algorithm usage without weights. The triplets are separated by spaces.

2. Alternative supertree. - Optional. Supertree of leaves that can be compared with the output tree. Should be in Newick format (see above).

3. Numerical parameters (Optional):
	-nd = maximal number of triplets to simulate (integer, >0, =10000000 by default). Please note that size of file with triplets should be <2Gb, hence number of triplets should be no more than 50000000.
	-nl = number of leafs in simulated binary tree (integer, >0)
	-mr = mutation rate for rooted triplets in simulation (value from 0 to 1, should be written using decimal dot ".", not comma ",", =0 by default)
	-w = usage of weights (0 for "Yes", 1 for "No"), =0 by default
	-index = index of program usage:
		1 - simulate tree and rooted triplets
		2 - make supertree based on rooted triplets
		3 - comparison of trees based on triplets and RF
		0 - 1+2+3

4. Names of files (Optional)
	-fit = file name with initial tree(s), no by default
	-fid = file name initial triplets, no by default
	-flg = file name log (="treeReconstruction.log" by default)
	-frt = file name with reconstructed tree (="resTree.tree" by default)
	-frtN = file name with reconstructed tree (numbers instead of names, ="res.tree" by default)
	-fsd = file name with triplets from simulated tree(s) (="qqq.dat" by default)
	-fst = file name with simulated tree (="simul.tree" by default)
	-fsg = file name with alternative tree, no by default

5. Output triplets - generated in the case of absence of input triplets. The format of this file is like for input triplets (see above).

6. Output supertree - supertree in Newick format (see above). The tree is rooted, non weighted and without distances.

Examples of input and output file are supplied with the tar file.
 
Running of software (note that in some computers you don't need to type "./" before the command):

1. Input the set of trees, generate corresponding triplets and construct supertree

Command: 
./treeFromTriplets -fit <fit> [-flg <flg>][-frt <frt>] [-frtN <frtN>][-fsd <fsd>][-nd <nd>]

Example:
./treeFromTriplets -fit example1.tree -nd 1000

For this example, resulted supertree can be found in file resTree.tree. Triplet-based and RF-based comparison with initial trees are summarised in file treeReconstruction.log.

2. Input the set of trees, generate corresponding triplets, construct supertree and compare with existing supertree

Command: 
./treeFromTriplets -fit <fit> [-flg <flg>][-frt <frt>] [-frtN <frtN>][-fsd <fsd>][-nd <nd>] -fsg <fsg>

Example:
./treeFromTriplets -fit example1.tree -nd 1000 -fsg example2.tree

For this example, resulted supertree can be found in file resTree.tree. Triplet-based and RF-based comparison with initial trees and with given supertrees are summarised in file treeReconstruction.log.

3. Calculate rooted triplet-based distance between two rooted trees (for given set of rooted triplets, no supertree construction):

Command: 
./treeFromTriplets -fit <fit> -fid <fid> [-flg <flg>] -fsg <fsg>[-w <w>] -index 3

Example:
./treeFromTriplets -fit example1.tree -fid qqq.dat -fsg example2.tree -index 3

4. Calculate rooted triplet-based distance between two rooted trees (no given set of rooted triplets, no supertree construction):

Command: 
./treeFromTriplets -fit <fit> -fid <fid> [-flg <flg>] -fsg <fsg>[-w <w>] -index 3

Example:
./treeFromTriplets -fit example1.tree -fid qqq.dat -fsg example2.tree -index 3

5. Simulations: Simulate random binary tree with nl leaves, randomly select up to nd rooted triplets with rate of mutation nr and reconstruct tree from the resulted list of triplets by using MaxCut-based algorithms 

Command: 
./treeFromTriplets -nl <nl> [-nd <nd>][-mr <mr>][-flg <flg>]

Example:
./treeFromTriplets -nl 50 -nd 100000 -mr 0.1

If you will meet any problem please contact with authors:
Zeev Frenkel, e-mail: zvfrenkel@gmail.com
Sagi Snir, e-mail: ssagi@research.haifa.ac.il
