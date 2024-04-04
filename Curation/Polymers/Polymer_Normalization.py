# this code is contains the functions needed to normalize polymers with the psmiles package

from psmiles import PolymerSmiles as PS


# the point of this is to create a wrapper for the psmiles package


class polymer_normalizer():
    def __init__(self) -> None:
        pass
    
    def normalize(psmiles):
        assert type(psmiles) == str, "PSMILES should be a string"
        if "*" not in psmiles:
            raise ValueError( "The enetered string appears not to be in PSMILES format")
        PS(psmiles).canonicalize 