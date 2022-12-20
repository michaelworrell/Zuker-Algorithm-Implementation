"""
Fills in recursive case entries in all tables.
"""
from TracebackObj import *
from freeEnergy import *


def fillWsubset(W, V, minLoop, ijDiff, sequence):
  """
  Fill in a slice of entries in table W where
  j - i equals ijDiff
  """
  n = len(sequence)

  # For each WM entry in subset, fill it in.
  for i in range(0, n - ijDiff):
    # Set j
    j = i + ijDiff

    # Case 1: j unpaired
    case1 = W[i][j-1][0]

    # Case 2: j paired
    case2 = float('inf')
    minimizerK = float('inf')

    # Dynamic programming solution to idenify the index along the sequence
    # which -- when paired with the nucleotide at index j -- results in
    # minimal free energy (given the minimum loop length constraint).
    for k in range(i, j - minLoop):
      partitionW = W[i][k - 1][0]
      if (i == k):
        partitionW = 0
      partitionV = V[k][j][0]
      if (partitionW + partitionV < case2):
        minimizerK = k
        case2 = partitionW + partitionV

    # Store information about minimum case
    W[i][j][0] = min(case1, case2)
    if case1 < case2:
      W[i][j][1].setCase(1)
    else:
      W[i][j][1].setCase(2)
      W[i][j][1].setK(minimizerK)

  # Return updated table
  return W


def fillVsubset(V, WM, ijDiff, sequence, energyTables):
  """
  Fill in a slice of entries in table W where
  j - i equals ijDiff. Uses energy tables to
  calculate free energy contributions of varior loop types.
  """
  n = len(sequence)

  # Fill in entry at V
  for i in range(0, n - ijDiff):
    j = ijDiff + i

    # Case 1: Hairpin Loop
    if checkComplementary(sequence[i],sequence[j]):
      eH = findEH(i, j, sequence, energyTables)
      case1 = eH
    else:
      case1 = float('inf')

    # Case 2: Stacking Loop
    if checkComplementary(sequence[i],sequence[j]):
      eS = findES(i, j, sequence, energyTables)
      case2 = V[i+1][j-1][0] + eS
    else:
      case2 = float('inf')

    # Calculate Case 3: Interior Loop / Bulge
    if checkComplementary(sequence[i],sequence[j]):
      minimizerI = float('-inf')
      minimizerJ = float('-inf')

      case3 = float('inf')
      for iPrime in range(i+1, j-1):
        for jPrime in range(iPrime + 1, j):
          eL = findEL(i, j, iPrime, jPrime, sequence, energyTables)
          if V[iPrime][jPrime][0] + eL < case3:
            case3 = V[iPrime][jPrime][0] + eL
            minimizerI, minimizerJ = iPrime, jPrime
    else:
      case3 = float('inf')

    # Caclulate Case 4: Multi-loop simplification
    if checkComplementary(sequence[i],sequence[j]):
      case4 = float('inf')
      minimizerK = float('inf')
      # Note: we did not implement a function that successfully calculated A for this
      # implementation of Zuker's Algorithm, but made it so that those functions
      # could be written and utilized in the future. For now, it returns a constant.
      a = getA(i, j, sequence)
      for k in range(i + 1, j):
        if (WM[i+1][k][0] + WM[k+1][j-1][0] + a < case4):
            case4 = WM[i+1][k][0] + WM[k+1][j-1][0] + a
            minimizerK = k
    else:
      case4 = float('inf')

    # Set WM[i][j] equal to the minimum of all
    # valid cases.
    caseMinimum = min(case1, case2, case3, case4)
    V[i][j][0] = caseMinimum

    # Store information about minimum case
    if caseMinimum == case1:
      V[i][j][1].setCase(1)
    elif caseMinimum == case2:
      V[i][j][1].setCase(2)
    elif caseMinimum == case3:
      V[i][j][1].setCase(3)
      V[i][j][1].setijprime(minimizerI, minimizerJ)
    else:
      V[i][j][1].setCase(4)
      V[i][j][1].setK(minimizerK)

  # Return updated table
  return V


def fillWMsubset(V, WM, ijDiff, sequence, energyTables):
  """
  Fill in a slice of entries in table WM where
  j - i equals ijDiff.
  """

  n = len(sequence)

  # Fill in entry at WM
  for i in range(0, n - ijDiff):
    # Get value of J
    j = ijDiff + i

    # Note: we did not implement functions that successfully calculated B or C
    # for this implementation of Zuker's Algorithm, but made it so that those functions
    # could be written and utilized in the future. For now, they return constants.
    b, c = getB(i, j, sequence), getC(i, j, sequence)

    # Case 1: j unpaired
    case1 = WM[i][j-1][0] + c
    # Case 2: i unpaired
    case2 = WM[i+1][j][0] + c
    # Case 3: Closed
    case3 = V[i][j][0] + b

    # Calculate Case 4: non-closed
    case4 = float('inf')
    minimizerK = float('inf')
    for k in range(i + 1, j):
      if (WM[i][k][0] + WM[k+1][j][0] < case4):
          minimizerK = k
          case4 = WM[i][k][0] + WM[k+1][j][0]

    # Set value at i, j to equal the minimum of all
    # valid cases
    caseMinimum = min(case1, case2, case3, case4)
    WM[i][j][0] = caseMinimum

    # Store information about minimum case
    if caseMinimum == case1:
      WM[i][j][1].setCase(1)
    elif caseMinimum == case2:
      WM[i][j][1].setCase(2)
    elif caseMinimum == case3:
      WM[i][j][1].setCase(3)
    else:
      WM[i][j][1].setCase(4)
      WM[i][j][1].setK(minimizerK)

  # Return updated table
  return WM
