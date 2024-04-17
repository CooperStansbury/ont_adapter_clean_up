import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.linalg import toeplitz
import cooler

import os
import sys


def dropZeroRows(A, t=0, return_ind=False):
    """A function to remove the zero count columns and rows from
    a symmetric matrix 
    
    args:
        : A (np.array): a symmetric matrix
        : t (int): number of contacts necesary to keep
        : return_ind (bool): if true, return the indices of dropped rows
    
    returns:
        : Ahat (np.array): a matrix with REDUCED dimensionality
    """
    rowSums = A.sum(axis=0)
    rmInd = np.argwhere(rowSums <= t)
    
    Ahat = A.copy()
    
    Ahat = np.delete(Ahat, rmInd, axis=0)
    Ahat = np.delete(Ahat, rmInd, axis=1)
    
    if return_ind:
        return Ahat, rmInd
    else:
        return Ahat
    
    
def forceAdjacentConnections(A, num=1):
    """A function to ensure that all adjcent genomic 
    loci are connected. Input is assumed to be binary. 
    
    args:
        : A (np.array): unweighted adjacency matrix
        : num (float): the value to set the off diagonal to
    
    returns:
        : Ahat (np.array): unweighted adjacency matrix with all i, i+1
        connections added 
    """
    Ahat = A.copy()
    
    for i in range(len(Ahat) - 1):
        Ahat[i, i+1] = num
        Ahat[i+1, i] = num
    return Ahat


def getToeplitz(A, return_means=False):
    """A function to get the toeplitz from the observed matrix 
    
    args:
        : A (np.array): the contact map
        : return_means (bool): if true, return the means
    
    returns:
        : E (np.array): the expected contact map
    """
    
    muDiags = []
    
    for offset in range(len(A)):
        # get each diagonal, divide it by it's
        # mean value and add it to the zero matrix
        mudiag = np.mean(np.diagonal(A, offset=offset))
        muDiags.append(mudiag)
        
    E = toeplitz(muDiags, muDiags)
    
    if return_means:
        return E, muDiags
    else:
        return E


def normalizeToeplitz(O):
    """A function to normalize and input matrix by the 
    mean diagonal value 
    
    args:
        : O (np.array): the 2d observed matrix to normalize
    
    returns:
        : A (np.array) of the same shape normalized by the diagonals
    """
    
    E = getToeplitz(O)
    A = np.divide(O, E)
    
    # handle NaNs
    A = np.where(np.isnan(A), 0, A)
    return A


def Abin(A, t=0):
    """A function to binarize a matrix by a threshold
    
    args:
        : A (np.array): matrix to binarize
        : t (int): a threshold
        
    returns:
        : Ahat (np.array): a binary matrix
    """
    return np.where(A > t, 1, 0)


def getIndOrder(clr, chromOrder):
    """A function to get indices of a resorted Hi-C mat
    args:
        : clr (cooler.Cooler): the cooler object
        : chromOrder (list): list of chromosome names in order
    returns:
        : newIndex (list of int): the new index order
    """
    newIndex = []    
    for chrom in chromOrder:
        start, end = clr.extent(chrom)
        subsetIndex = list(range(start, end))
        newIndex += subsetIndex
    return newIndex
    

def reOrgCooler(clr, newIndex, balance=True):
    """A function to index a matrix based on a chromosome 
    order
    
    args:
        : clr (cooler.Cooler): the cooler object
        : newIndex (list of int): the new index order
        : balance (bool): if true, balance the matrix
    returns:
        : HiC (np.array): reindexed by the chromosome order
    """
    mat = clr.matrix(balance=balance)[:]
    HiC = mat[:, newIndex][newIndex, :]
    return HiC