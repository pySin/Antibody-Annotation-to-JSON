# Populate JSON file with antibody annotation data
import json


class AntibodyToJSON:

    def __init__(self, path):
        self.path = path


def main():
    antibody_json = AntibodyToJSON("samples")
