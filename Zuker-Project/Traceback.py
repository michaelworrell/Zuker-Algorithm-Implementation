"""
Constructs the optimal set of pairings using traceback
"""
from Traceback import *


def tracebackWM(V, WM, i, j):
  """
  Performs traceback based on cases in table WM
  """
  tracebackObj = WM[i][j][1]
  if tracebackObj.case == 1:
    # Performs traceback given that Case 1 was optimal
    pairings = tracebackWM(V, WM, i, j-1) + ["."]
  elif tracebackObj.case == 2:
    # Performs traceback given that Case 2 was optimal
    pairings = ["."] + tracebackWM(V, WM, i+1, j)
  elif tracebackObj.case == 3:
    # Performs traceback given that Case 3 was optimal
    pairings = tracebackV(V, WM, i, j)
  else:
    # Performs traceback given that Case 4 was optimal
    k = tracebackObj.k
    pairings = tracebackWM(V, WM, i, k) + tracebackWM(V, WM, k+1, j)

  # return partial pairings list.
  return pairings


def tracebackV(V, WM, i, j):
  """
  Performs traceback based on cases in table V
  """
  tracebackObj = V[i][j][1]
  
  if tracebackObj.case == 1:
    # Tracback case 1
    length = j - i + 1
    pairings = ["." for i in range(length)]
    pairings[0] = "("
    pairings[length - 1] = ")"
  elif tracebackObj.case == 2:
    # Tracback case 2
    pairings = tracebackV(V, WM, i+1, j-1)
    pairings = ["("] + pairings + [")"]
  elif tracebackObj.case == 3:
    # Tracback case 3
    iprime = tracebackObj.iprime
    jprime = tracebackObj.jprime
    # Get right, middle, and left sides
    left = ["." for i in range(iprime- i)]
    middle = tracebackV(V, WM, iprime, jprime)
    right = ["." for i in range(j- jprime)]
    # Add pairing params
    left[0] = "("
    right[j- jprime - 1] = ")"
    # Merge together
    pairings = left + middle + right
  elif tracebackObj.case == 4:
    # Tracback case 4
    k = tracebackObj.k
    partition1 = tracebackWM(V, WM, i+1, k)
    partition2 = tracebackWM(V, WM, k+1, j-1)
    pairings = ["("] + partition1 + partition2 + [")"]

  # return partial pairings list.
  return pairings


def tracebackW(W, V, WM, i, j):
  """
  Performs traceback based on cases in table W
  """
  # Handle Edge Case
  if i > j:
    return []
  
  tracebackObj = W[i][j][1]

  if tracebackObj.case == 1:
    # Performs traceback given that Case 1 was optimal
    partition1 = tracebackW(W, V, WM, i, j-1)
    pairings = partition1 + ["."]
  elif tracebackObj.case == 2:
    # Performs traceback given that Case 2 was optimal
    k = tracebackObj.k
    partition1 = tracebackW(W, V, WM, i, k-1)
    partition2 = tracebackV(V, WM, k, j)
    pairings = partition1 + partition2
  else:
    pairings = ["." for i in range(j - i + 1)]
  
  # return partial pairings list.
  return pairings
