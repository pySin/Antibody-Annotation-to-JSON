# Populate JSON file with antibody annotation data
import json
import os


class AntibodyToJSON:

    def __init__(self, path):
        self.path = path
        self.files = os.listdir(path)

    def read_define_records(self, filename):
        with open(filename, "r") as f:
            content = [record.strip() for record in f.read().split(";")]
            # for item in content:
            #     print(item)
            print(content)


def main():
    antibody_json = AntibodyToJSON("samples")
    antibody_json.read_define_records(antibody_json.files[0])


if __name__ == "__main__":
    main()
