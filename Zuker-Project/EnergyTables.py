"""
Class that stores free energy scoring tables and lists, 
and uses them to determine free energy contributions
for various types of loops.
"""
import sys

class EnergyTables:
  """
  Initialize the tables used to calculate free energy.
  """
  def __init__(self, sequence):
    self.stacking_Table = energy_table_from_file("stackingTable.txt", "Stacking")
    self.sequence = sequence
    self.internalLength, self.bulgingLength, self.hairpinLength = loop_tables_from_file("loopLength.txt")
    self.hairpin_Table = energy_table_from_file("terminalMismatchTable.txt", "TerminalMismatch")

  """
  Get energy contribution from stacking loop (eL)
  """
  def get_stacking_energy(self, bottom, top):
    
    i = self.stacking_Table[0].index(bottom.upper())
    j = self.stacking_Table[0].index(top.upper())
    energy_contribution = self.stacking_Table[i][j]
    # print("Free Energy = {0}".format(energy_contribution))
    return energy_contribution

  """
  Get energy contribution from hairpin loop (eH)
  """
  def get_hairpin_energy(self, string):
    energy_contribution = self.hairpinLength[len(string) - 2];
    
    bottom = string[0] + string[-1]
    top = string[1] + string[-2]
    i = self.hairpin_Table[0].index(bottom.upper())
    j = self.hairpin_Table[0].index(top.upper())

    energy_contribution += self.hairpin_Table[i][j];
    return energy_contribution
  

def energy_table_from_file(fileName, tableToIsolate):
  """
  Read an energy table from a file, assuming similar
  loop scoring tables share an idential format except
  for the label preceeding the matrix information.
  """
  with open(fileName, "r") as fileToRead:
    s = fileToRead.read()
    s = s.split("\n")
    fileToRead.close()

    # Read through lines until table is about to start.
    tableLabelEncountered = False
    while tableLabelEncountered == False:
      if s[0] == tableToIsolate:
        tableLabelEncountered = True
      s = s[1:]

    # Initialize the scoring matrix
    label_list = [" "] + s[0].split(" ")
    stack_table = ["" for i in range(len(label_list))]
    stack_table[0] = label_list
    s = s[1:]

    # Fill in the scoring matrix from the file
    for i in range(1, len(stack_table[0])):
      stack_table[i] = s[0].split(' ')
      for j in range(1, len(stack_table[0])):        
        if stack_table[i][j] == ".":
          stack_table[i][j] = float('inf')
        else:
          stack_table[i][j] = float(stack_table[i][j])
      s = s[1:]
  
  # return the completed scoring matrix table
  return stack_table
 

def loop_tables_from_file(fileName):
  """
  Read in the loop length initialization score
  for internal, bulging, and hairpin loops and
  returns the information in the form of lists. 
  """
  with open(fileName, "r") as fileToRead:
    s = fileToRead.read()
    s = s.split("\n")
    fileToRead.close()
    header = 'LoopLength30'
    headerEncountered = False
    
    # Find head of the tables
    while headerEncountered == False:
      if s[0] == header:
        headerEncountered = True
      s = s[1:]

    # internalize lists
    internal_list = ['internal']
    bulging_list = ['bulging']
    hairpin_list = ['hairpin']
    s = s[1:]

    # Adds entries from file to appropriate lists.
    for i in range(30):
      entry = s[0].split(' ')
      # Add entry to internal_loop_list
      if entry[0] == ".":
        internal_list.append(float('inf'))
      else:
        internal_list.append(float(entry[0]))
      # Add entry to bulge_loop_list
      if entry[1] == ".":
        bulging_list.append(float('inf'))
      else:
        bulging_list.append(float(entry[1]))
      # Add entry to hairpin_loop_list
      if entry[2] == ".":
        hairpin_list.append(float('inf'))
      else:
        hairpin_list.append(float(entry[2]))
      s = s[1:]
  
  # return the lists.
  return internal_list, bulging_list, hairpin_list
