import sys
from TracebackObj import *
from fillSubsets import *
from traceback import *
from EnergyTables import *

def read_single_fasta(filename):
    """
    Reads in a single fasta file (only contains one sequence) and returns
    the sequence from the file.
    Input: filename (String) - the name of the fasta file to read.  The
            sequence may be on multiple lines in the fasta file.
    Output: sequence (String) - a single string representing the sequence in
            the fasta file.
    """
    # Scans the file until appropriate section is reached.
    with open(filename) as fileToRead:
    	firstLineEncountered = False
    	
    	while not firstLineEncountered:
    		lineToRead = fileToRead.readline()
    		if (lineToRead[0] == ">"):
    			firstLineEncountered = True
    	
    	# Read the sequence
    	single_fasta = fileToRead.read()
    	
    # Removes line feed from the code
    single_fasta = single_fasta.replace('\n', '')
    # Converts the string to Uppercase
    single_fasta = single_fasta.upper()
    
    return single_fasta


def initializeTables(sequence, minloop):
  """
  Initialize all of the tables (W, V, and WM) by filling in
  the base cases, and filling in -inf for entries that will be
  filled out recursively 
  """
  length = len(sequence)

  # - infinity initialization for bugfixing purposes
  W = [[[float('-inf'), Traceback()] for i in range(length)] for j in range(length)]
  V = [[[float('-inf'), Traceback()] for i in range(length)] for j in range(length)]
  WM = [[[float('-inf'), Traceback()] for i in range(length)] for j in range(length)]
  
  # Three n*n matrices are needed. W[i][j] stores information about whether the ith base in the sequence is paired to the jth base in the sequence. V[i][j] stores information about what kind of structures are built from ith base and jth base. WM[i][j] stores information about a particular kind of structure--the inner loop--because it is complicated. All three matrices would be used for traceback.
  for i in range(length):
    for j in range(length):
      # Loops that are too small are energetically unfavorable. Therefore, a limit on loop length helps us initialize the table.
      if j <= i + minloop:
        W[i][j][0] = 0
        V[i][j][0], WM[i][j][0] = float('inf'), float('inf')
  return W, V, WM



def fillTables(W, V, WM, minLoop, sequence):
  """
  Fill out the recursive case entries for each table
  """
  n = len(sequence)
  # Get the tables containing free energy scoring list and matricies.
  energyTables = EnergyTables(sequence)

  # Main Fill Loop for V and WM
  # Filling in slices of entries in V and WM
  for ijDiff in range(minLoop + 1, n):
      # Fill in entries in V such that j - i = ijDiff
      V = fillVsubset(V, WM, ijDiff, sequence, energyTables)
      # Fill in entries in WM where j - i = ijDifff 
      WM = fillWMsubset(V, WM, ijDiff, sequence, energyTables)
  
  # Main Fill Loop for W
  for ijDiff in range(minLoop + 1, n):
      W = fillWsubset(W, V, minLoop, ijDiff, sequence)

  # return filled tables
  return W, V, WM


def traceback(W, V, WM, sequence):
  """
  Perform traceback on the tables and return a parenthases and dots
  representation of the pairs present in the folding.
  """
  n = len(sequence)
  return tracebackW(W, V, WM, 0, n - 1);



def main():
  """
  Main Function Performs Zuker's Algorithm
  """
  # Gets parameters from command line
  args = sys.argv[1:]
  if len(args) != 2:
    print("Improper number of arguments provided.")
    return
  fileName = args[0]
  minLoop = int(args[1])

  # Read sequence from file
  sequence = read_single_fasta(fileName)
  # Convert T's to U's (if DNA sequence provided instead of RNA)
  sequence = sequence.replace("T", "U")
  
  # Fills out the tables in their entirety
  W, V, WM =initializeTables(sequence,minLoop)
  W, V, WM =fillTables(W, V, WM, minLoop, sequence)
  
  # Performs tracebakc to get a string representation of pairings
  optimalFolding = tracebackW(W, V, WM, 0, len(sequence) - 1)
  stringFold = "".join(optimalFolding)

  # Prints the result
  print("Total Free Energy: {0}".format(W[0][len(sequence) - 1][0]))
  print(stringFold)

  
if __name__ == "__main__": 
  main()
