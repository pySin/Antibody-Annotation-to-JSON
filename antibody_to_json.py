#  Populate JSON file from text data


class AntibodyTxtJSON:

    def __init__(self):
        self.test = None

    def txt_to_dict(self, filename):
        with open(filename, "r") as f:
            for line in f.read():
                print(line)


def main():
    AntibodyJSON = AntibodyTxtJSON
    AntibodyJSON.txt_to_dict()
