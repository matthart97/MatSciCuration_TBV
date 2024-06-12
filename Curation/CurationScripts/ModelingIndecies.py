import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.typing import ArrayLike, NDArray
import rogi


"""

TODO this is code written by the cursed language model, fix it @MATT

"""



def calculate_modi(distance_matrix: NDArray, target: ArrayLike) -> float:
    """
    Calculates the Modelability Index (MODI) of a dataset.

    Args:
        distance_matrix: Numpy array of the distance matrix in square form.
        target: Target array. N.B. Target must be binary!

    Returns:
        The MODI value.
    """
    try:
        modi_index_value = rogi.MODI(Dx=distance_matrix, Y=target)
        return modi_index_value
    except Exception as e:
        print(f'There was an error with your input: {e}')
        return None

def calculate_rogi(target: ArrayLike, smiles: ArrayLike, use_descriptors=True, descriptors: pd.DataFrame = None) -> float:
    """
    Calculates the Roughness Index (ROGI) of a dataset.

    Args:
        target: Column of targets from dataframe.
        smiles: Column of smiles from dataframe.
        use_descriptors: Option to use precomputed descriptors. Default is False, in this case
                         Morgan fingerprints of length 2048 & radius 2 are computed. If True, user
                         passes the dataframe of descriptors excluding the target column.
        descriptors: Valid pandas dataframe containing calculated descriptors.
                     default = None

    Returns:
        The ROGI value.
    """
    if use_descriptors:
        assert isinstance(descriptors, pd.DataFrame), "Descriptors must be a pandas DataFrame when use_descriptors is True."

    try:
        if use_descriptors:
            ri = rogi.RoughnessIndex(Y=target, X=descriptors, metric='euclidean')
        else:
            ri = rogi.RoughnessIndex(Y=target, smiles=smiles)
        rogi_index_value = ri.compute_index()
        return rogi_index_value
    except Exception as e:
        print(f'There was an error with your input: {e}')
        return None

def calculate_sari(distance_matrix: NDArray, target: ArrayLike) -> float:
    """
    Calculates the Structure-Activity Relationship Index (SARI) of a dataset.

    Args:
        distance_matrix: Numpy array of the distance matrix in square form.
        target: Target array.

    Returns:
        The SARI value.
    """
    try:
        # Example of a potential SARI calculation - this should be replaced with the actual method
        # As SARI is a hypothetical index here, let's assume we have a function in rogi called SARI
        sari_index_value = rogi.SARI(Dx=distance_matrix, Y=target)
        return sari_index_value
    except Exception as e:
        print(f'There was an error with your input: {e}')
        return None

def calculate_square_distance_matrix(descriptors: pd.DataFrame = None) -> NDArray:
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
            print(f'There was an error with your input: {e}')
            return None
    else:
        print("Descriptor dataframe not provided.")
        return None


# Example usage:

# Load the descriptor data
df = pd.read_csv('/path/to/DyeDescriptors.csv')

# Get the descriptors and target properties
descriptors = df.filter(regex='fp')
target_property = df['Absorption max (nm)'].values  # Example target property
smiles_column = df['SMILES']

# Calculate the distance matrix
distance_matrix = calculate_square_distance_matrix(descriptors)

# Calculate ROGI
rogi_value = calculate_rogi(target=target_property, smiles=smiles_column, descriptors=descriptors)
print(f'ROGI value: {rogi_value}')

# Calculate MODI (assuming the target property is binary for this example)
modi_value = calculate_modi(distance_matrix=distance_matrix, target=target_property)
print(f'MODI value: {modi_value}')

# Calculate SARI
sari_value = calculate_sari(distance_matrix=distance_matrix, target=target_property)
print(f'SARI value: {sari_value}')
