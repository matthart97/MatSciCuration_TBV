import numpy as np
from sklearn.ensemble import IsolationForest
from sklearn.svm import OneClassSVM
from sklearn.neighbors import LocalOutlierFactor

class NoveltyDetector:
    def __init__(self, method='iforest', random_state=None):
        self.method = method
        self.random_state = random_state
        self.model = self._select_model()

    def _select_model(self):
        if self.method == 'iforest':
            return IsolationForest(random_state=self.random_state)
        elif self.method == 'ocsvm':
            return OneClassSVM()
        elif self.method == 'lof':
            return LocalOutlierFactor(novelty=True)
        else:
            raise ValueError(f"Unknown method: {self.method}")

    def fit(self, X):
        if self.method in ['iforest', 'ocsvm']:
            self.model.fit(X)
        else:
            raise ValueError("LOF does not support fit method for novelty detection. Use fit_predict.")

    def predict(self, X):
        if self.method == 'lof':
            raise ValueError("LOF cannot predict on new data directly. Use fit_predict for training data.")
        return self.model.predict(X)

    def fit_predict(self, X):
        if self.method == 'lof':
            return self.model.fit_predict(X)
        else:
            self.fit(X)
            return self.predict(X)