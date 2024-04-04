
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.typing import ArrayLike, NDArray


import rogi


def rogi_index(target, smiles,descriptors: pd.DataFrame = None, use_descriptors: bool = False):
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


df = pd.read_csv("/home/matt/Proj/MatSciCuration_TBV/Data/Polymers/PolymerDescriptors.csv")

Desc = df.filter(regex = 'fp')

# make a list of properties that we need to calulate the rogi index for

