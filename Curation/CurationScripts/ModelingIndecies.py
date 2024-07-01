import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.typing import ArrayLike, NDArray
import rogi


"""

TODO this is code written by the cursed language model, fix it @MATT

"""



def nearest_neighbors(reference, query, k=1, self_query=False, return_distance=False):
    """
    Gets the k nearest neighbors of reference set for each row of the query set

    Parameters
    ----------
        reference : array_like
            An array of points to where nearest neighbors are pulled.
        query : array_like
            An array of points to query nearest neighbors for
        k : int > 0, optional
            the number of nearest neighbors to return
        self_query : bool, optional
            if reference and query are same set of points, set to True
            to avoid each query returning itself as its own nearest neighbor
        return_distance : bool, optional
            if True, return distances of nearest neighbors
    Returns
    -------
        i : integer or array of integers
            The index of each neighbor in reference
            has shape [q, k] where q is number of rows in query
        d : float or array of floats, optional
            if return_distance set to true, returns associated
            euclidean distance of each nearest neighbor
            d is element matched to i, ei the distance of i[a,b] is d[a,b]
            The distances to the nearest neighbors.
    """

    tree = sp.KDTree(reference)

    if self_query:
        k = [x+2 for x in range(k)]
    else:
        k = [x+1 for x in range(k)]

    d, i = tree.query(query, k=k, workers=-1)

    if return_distance:
        return i, d
    else:
        return i


def modi(data, labels, return_class_contribution=False):
    """
    Gets the MODI from the given data and label set

    Parameters
    ----------
        data : array_like
            An array chemical descriptors (rows are chemicals and columns are descriptors).
        labels : array_like
            An array labels that are row matched to the data array
        return_class_contribution : bool, optional
            if True, return the normalized MODI for each class. Useful for imbalanced datasets
    Returns
    -------
        modi : float
            the calculated MODI for the given data and label
        class_contrib : list of tuples of length 2 (str, float), optional
            if return_class_contribution set to true, returns associated
            MODI score for each class in the data as a tuple of (class, MODI)
    """
    # get all the classes present in the dataset
    classes = np.unique(labels)
    k = classes.shape[0]

    # get the labels of the nearest neighbors
    nn_idx = nearest_neighbors(data, data, k=1, self_query=True)
    nn_labels = labels[nn_idx]

    # calculate the modi
    modi_value = 0
    class_contrib = []

    # loop through each class
    for c in classes:
        c_arr = np.where(labels == c)[0]
        c_labels = labels[c_arr]
        c_nn_labels = nn_labels[c_arr].flatten()

        modi_value += np.sum(c_labels == c_nn_labels) / c_arr.shape[0]
        class_contrib.append((c, np.sum(c_labels == c_nn_labels) / c_arr.shape[0]))

    if not return_class_contribution:
        return (k ** -1) * modi_value
    else:
        return (k ** -1) * modi_value, class_contrib


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
