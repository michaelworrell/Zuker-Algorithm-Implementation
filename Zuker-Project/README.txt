Different files were built to handle different components of Zuker algorithm. More details can be found in these files separately.

To run Zuker, on the command line, call the file, input the name of the fasta file, and user-defined minmal loop length.

1. Filling Tables

This is handled by fillSubsets.py. It contains code that fills out free energy in slices in table W, V and WM as well as initializing traceback objects also stored in the tables.

2.Energy Contribution

This is handled by EnergyTables.py and freeEnergy.py. EnergyTables.py stores free energy scoring tables and lists, and uses them to determine free energy contributions for various types of loops. There are tables for hairpin, stacking and internal loop and only constants for multiloops. This is because energy calculation for multiloops is simplified in Zuker, or the algorithm would take on exponential run time. The user can input different parameters for energy calculation by replacing our table with user-defined ones of the same format. Our energy tables (stackingTable.txt, stackingTable(copy)1.txt, terminalMismatchTable.txt, loopLength.txt) are adapted from Mathews, D.H., Disney, M.D., Childs, J.L., Schroeder, S.J., Zuker, M. and Turner, D.H. (2004) Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proc. Natl. Acad. Sci. USA, 101, 7287-7292. We obtained information from the website
https://rna.urmc.rochester.edu/NNDB/turner04/index.html

freeEnergy.py calculates the energy of an input loop given the stored energy information. 

3. Traceback

This is handled by traceback.py and Traceback.py. Traceback.py defines the Traceback objects we used to store traceback information. traceback.py contains code that retrieve the optimal structure computed from filled W, V and WM tables.

4. RNA secondary structure prediction

zukers.py wraps all the components together. It initializes the tables, fills the tables slice by slice, and performs traceback. To run Zuker, on the command line, call the file, input the name of the fasta file, and user-defined minmal loop length. The T bases in the input fasta file would be converted to U for predicting the structure of the RNA of the specified sequence.

5.sampleData.fasta

Here we included the sequence of MT-TI as the sample data.
