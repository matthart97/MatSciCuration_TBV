import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from rogi import SARI, MODI

class ModelabilityChecker:
    def __init__(self, data, labels):
        self.data = data
        self.labels = labels

    def knn_modelability(self, n_neighbors=3):
        model = KNeighborsClassifier(n_neighbors=n_neighbors)
        model.fit(self.data, self.labels)
        predictions = model.predict(self.data)
        accuracy = accuracy_score(self.labels, predictions)
        return accuracy

    def random_forest_modelability(self, n_estimators=100):
        model = RandomForestClassifier(n_estimators=n_estimators)
        model.fit(self.data, self.labels)
        predictions = model.predict(self.data)
        accuracy = accuracy_score(self.labels, predictions)
        return accuracy

    def sari_index(self, smiles):
        sari = SARI(pKi=self.labels, smiles=smiles)
        return sari.compute_sari()

    def modi_index(self, distance_matrix):
        return MODI(Dx=distance_matrix, Y=self.labels)

    def rogi_index(self, Y, smiles=None, descriptors=None, metric='euclidean'):
        from rogi import RoughnessIndex
        ri = RoughnessIndex(Y=Y, smiles=smiles, X=descriptors, metric=metric)
        return ri.compute_index()

# Usage example
# data = np.array([[...]])  # your dataset
# labels = np.array([...])  # your labels
# checker = ModelabilityChecker(data, labels)
# print(checker.knn_modelability())
# print(checker.random_forest_modelability())
# print(checker.sari_index(smiles_data))
# print(checker.modi_index(distance_matrix))
# print(checker.rogi_index(labels, smiles_data))
