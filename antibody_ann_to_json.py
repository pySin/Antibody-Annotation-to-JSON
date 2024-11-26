# Populate JSON file with antibody annotation data
import json
import os


class AntibodyToJSON:

    def __init__(self, path):
        self.path = path
        self.files = os.listdir(path)
        self.antibody_ann_dict = {}
        self.methods = {
            "Antigen": self.antigen_record
        }

    def read_devide_records(self, filename):
        with open(filename, "r") as f:
            # Split by semicolon and filter new lines(\n)
            content = [record.replace("\n", "") for record in f.read().split(";")]
            records = []
            # Split records further if they have double slash(//)
            for r in content:
                if "//" in r:
                    chains = r.split("//")
                    [records.append(c) for c in chains]
                else:
                    records.append(r)
            # [print(re) for re in records]
            return records

    def single_file_transfer(self, filename):
        # Produce JSON file from .txt annotation file
        current_records = self.read_devide_records(filename)
        for record in current_records:

            key = record.split(":")[0]
            key = key[:key.index("[")] if "[" in key else key

            if key in self.methods:
                self.methods[key](record)

    def antigen_record(self, record):
        key, value = record.split(":")
        self.antibody_ann_dict[key] = " ".join(value.split(" ")[:-1]).strip()
        self.antibody_ann_dict[key + "-Gene"] = value.split(" ")[-1][1:-1]
        print(f"Main Dictionary: {self.antibody_ann_dict}")

    def note_record(self, record):
        pass
