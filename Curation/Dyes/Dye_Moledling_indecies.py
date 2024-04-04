import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.typing import ArrayLike, NDArray


import rogi


def rogi_index( target, smiles, use_descriptors=True,descriptors= None,) -> float:
    """Calculates the roughness index (ROGI) of a dataset.

    Args:
      descriptors: Valid pandas dataframe containing calculated descriptors.
                   default = None
      target:  Column of targets from dataframe.
      smiles: Column of smiles from dataframe.
      use_descriptors: Option to use precomputed descriptors. Default is False, in this case
                       Morgan fingerprints of lenght 2048 & radius 2 are computed. If True, user
                       passes the dataframe of descriptors exluding the target column.

    Returns:
      The ROGI value.
    """
    assert type(descriptors) == pd.DataFrame
    try:
        if use_descriptors == True:
            ri = rogi.RoughnessIndex(Y=target, X=descriptors, metric='euclidean')
            rogi_index = ri.compute_index()
        else:
            ri = rogi.RoughnessIndex(Y=target, smiles=smiles)
            rogi_index = ri.compute_index()
        return rogi_index
    except Exception as e:
        print('there was an error with your input :{0}'.format(e))


def modi_index(distance_matrix: NDArray, target: ArrayLike) -> float:
    """Calculates the modelability index (MODI) of a dataset from a square distance maxtrix.
       N.B. Target must be binary! (https://pubs.acs.org/doi/10.1021/ci400572x)

    Args:
      distance_matrix: array of distance matrix of X in square form.
      target: Target array. N.B. Target must be binary!

    Returns:
      The MODI value.
    """
    try:
        modi_index = rogi.MODI(Dx=distance_matrix, Y=target)
        return modi_index
    except Exception as e:
        print('there was an error with your input :{0}'.format(e))


def rmodi_index(distance_matrix: NDArray, target: ArrayLike, delta: float = 0.625) -> float:
    """Calculates the regression modelability index (RMODI) of a dataset.
       (https://pubs.acs.org/doi/10.1021/acs.jcim.8b00313)

    Args:
      distance_matrix: array of distance matrix of X in square form.
      target: Target array.
      delta: tunable parameter

    Returns:
      The RMODI value.
    """
    try:
        rmodi_index = rogi.RMODI(Dx=distance_matrix, Y=target, delta=delta)
        return rmodi_index
    except Exception as e:
        print('there was an error with your input :{0}'.format(e))



def calculate_sqaure_distance_matrix(descriptors: pd.DataFrame = None) -> NDArray:
    """Calculates Pairwise distances between observations in descriptor dataframe

    Args:
      descriptors: Valid pandas dataframe containing calculated descriptors.
                   default = None

    Returns:
      The distance matrix in square form
    """
    if descriptors is not None:
        try:
            matrix = pdist(descriptors)
            square_matrix = squareform(matrix)
            return square_matrix
        except Exception as e:
            print('there was an error with your input :{0}'.format(e))
    else:
        print("Descriptor dataframe not provided.")



# load in the descriptor data
df = pd.read_csv('/home/matt/Proj/MatSciCuration_TBV/Data/Dyes/DyeDescriptors.csv')

Descriptors = df.filter(regex='fp')


# get the list of properties that we want

props=['Absorption max (nm)', 'Emission max (nm)',
       'Lifetime (ns)', 'Quantum yield', 'log(e/mol-1 dm3 cm-1)',
       'abs FWHM (cm-1)', 'emi FWHM (cm-1)', 'abs FWHM (nm)', 'emi FWHM (nm)',
       'Molecular weight (g mol-1)', 'vie', 'aie',
       'vea', 'aea', 'hl', 's0s1', 's0t1', 'hr', 'cr2', 'cr1', 'er', 'ar1',
       'ar2', 'lumo', 'homo']

#test the code with 


#rogi_index(df[props[0]],df['SMILES'],descriptors=Descriptors)



distance_matrix = calculate_sqaure_distance_matrix(Descriptors)
try:
  test = rmodi_index(distance_matrix,df[props[0]])
except: 
  print("what the FUCK")

with open("output.txt") as f:
    f.write(test)


# set up main loop 