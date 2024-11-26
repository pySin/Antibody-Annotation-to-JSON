# Populate JSON file with antibody annotation data
import json
import os


class AntibodyToJSON:

    def __init__(self, path):
        self.path = path
        self.files = os.listdir(path)
        # print(self.files)
        self.antibody_ann_dict = {}
        self.old_key = None
        self.methods = {
            "Antigen": self.antigen_record,
            "Note": self.note_record
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

            if record.startswith("Note"):
                self.note_record(record)
                continue

            key = record.split(":")[0]
            key = key[:key.index("[")] if "[" in key else key

            if key in self.methods:
                self.methods[key](record)
            self.old_key = key

    def antigen_record(self, record):
        key, value = record.split(":", 1)
        self.antibody_ann_dict[key] = " ".join(value.split(" ")[:-1]).strip()
        self.antibody_ann_dict[key + "-Gene"] = value.split(" ")[-1][1:-1]
        # print(f"Main Dictionary: {self.antibody_ann_dict}")

    def note_record(self, record):
        key, value = record.split(":", 1)

        if key != "Note":
            note_range = key.replace("Note", "")
            # print(f"Note range: {note_range}")
            self.antibody_ann_dict[self.old_key + note_range + "-Note"] = value.strip()
        else:
            self.antibody_ann_dict[self.old_key + "-" + key] = value.strip()
        for item in self.antibody_ann_dict:
            print(f"Antibody Dict Item: {item}")

