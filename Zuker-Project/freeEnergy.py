"""
Functions used to calculate free energy
"""
from EnergyTables import *


def checkComplementary(a,b):
  """
  Checks whether two bases can be paired to each other (includes) complementary and GU pairs. Returns true if they can be paired, false otherwise.
  """
  # Checks check whether two bases are pairable or not
  bases = [('A','U'),('C','G'),('G', 'U')]
  if ((a,b) in bases) or ((b, a) in bases):
    return True
  else:
    return False


def findEH(i, j, sequence, energyTables):
  """
  Find Energy Contribution of Hairpin Loop
  """
  string = sequence[i:j+1]
  # Assumes no hairpin loop with length greater than 30
  # (limitation of our scoring tables)
  if (len(string) - 2) > 30:
    return float('inf')
  # Find and return free energy contribution of hairpin loop.
  free_energy = energyTables.get_hairpin_energy(string)
  return free_energy

def findEL(i, j, iPrime, jPrime, sequence, energyTables):
  """
  Find Energy Contribution of Interior / Bulge Loop

  NOTE: Until Proper Methodoly implemented, we will assume that there
  are no internal or bulging loops in an optimal solution
  """
  return float('inf')

def findES(i, j, sequence, energyTables):
  """
  Find Energy Contribution of Stacking Loop
  """
  top = sequence[i + 1] + sequence[j - 1]
  bottom = sequence[i] + sequence[j]
  free_energy = energyTables.get_stacking_energy(bottom, top)
  return free_energy

def getA(i, j, sequence):
  # Get Energy Contribution for Closing of Loop (Weight)
  return 9.25
def getB(i, j, sequence):
  # Get Energy Contribution weight b
  return 0.91
def getC(i, j, sequence):
  # Get Energy Contribution weight c
  return -0.63

