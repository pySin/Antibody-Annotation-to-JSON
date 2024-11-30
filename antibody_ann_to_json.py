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

        for t_key, t_value in self.antibody_ann_dict.items():
            print(f"Key: {t_key}, Value: {t_value}")
        # print(self.antibody_ann_dict)

    def antigen_record(self, record):
        key, value = record.split(":", 1)

        if "[" in key:
            name_key, region = key.split("[")
            region = region[:-1]
            name = value.split(",")[0]
            gene = value.split(" ")[-1][1:-1]
            if "Antigen" not in self.antibody_ann_dict:
                antigen_data = [{"Region": region, "Name": name, "Gene": gene}]
                self.antibody_ann_dict[name_key] = antigen_data
            else:
                antigen_data = {"Region": region, "Name": name, "Gene": gene}
                self.antibody_ann_dict[name_key].append(antigen_data)

            # region = [int(num) for num in region[:-1].split(",")]
            # print(f"Name: {name}, Region: {region}")
            # tuple_key = (name, tuple(region))
            # self.antibody_ann_dict[tuple_key] = " ".join(value.split(" ")[:-1]).strip()

        else:
            self.antibody_ann_dict[key] = " ".join(value.split(" ")[:-1]).strip()

        # self.antibody_ann_dict[key + "-Gene"] = value.split(" ")[-1][1:-1]

    def note_record(self, record):
        key, value = record.split(":", 1)

        if key != "Note":
            note_range = key[5:-1]
            # print(f"Note range: {note_range}")
            # self.antibody_ann_dict[self.old_key + note_range + "-Note"] = value.strip()
        else:
            self.antibody_ann_dict[self.old_key + "-" + key] = value.strip()
        # for item in self.antibody_ann_dict:
        #     print(f"Antibody Dict Item: {item}")

