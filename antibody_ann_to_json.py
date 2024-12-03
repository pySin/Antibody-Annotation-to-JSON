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
            "Note": self.note_record,
            "Antigen-Note": self.antigen_note,
            "CDR": self.cdr_record,
            "Heavy Chain": self.heavy_chain_record,
            "Light Chain": self.light_chain_record,
            "Chain": self.chain_record
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
            records = [r for r in records if len(r) > 2]
            return records

    def single_file_transfer(self, filename):

        # Produce JSON file from .txt annotation file
        current_records = self.read_devide_records(filename)
        for record in current_records:

            if record.startswith("Note"):
                self.methods["Note"](record)
                continue

            if record.startswith("CDR"):
                self.methods["CDR"](record)
                continue

            key = record.split(":")[0]
            key = key[:key.index("[")] if "[" in key else key

            if key in self.methods:
                self.methods[key](record)
            else:
                self.normal_record(record)
            self.old_key = key

        # print(f"Full Dict: {self.antibody_ann_dict}")

        filename = filename.split("/")[1]
        with open(f"{filename.split('.')[0]}.json", "w") as jf:
            json.dump(self.antibody_ann_dict, jf, indent=4)

    def antigen_record(self, record):
        key, value = record.split(":", 1)

        if "[" in key:
            name_key, instance = key.split("[")
            instance = instance[:-1].strip()
            name = value.split(",")[0].strip()
            gene = value.split(" ")[-1][1:-1].strip()
            if "Antigen" not in self.antibody_ann_dict:
                antigen_data = [{"Instance": instance, "Name": name, "Gene": gene}]
                self.antibody_ann_dict[name_key] = antigen_data
            else:
                antigen_data = {"Instance": instance, "Name": name, "Gene": gene}
                self.antibody_ann_dict["Antigen"].append(antigen_data)

        else:
            self.antibody_ann_dict[key] = " ".join(value.split(" ")[:-1]).strip()

        # self.antibody_ann_dict[key + "-Gene"] = value.split(" ")[-1][1:-1]

    def note_record(self, record):
        key, value = record.split(":", 1)

        if key != "Note":
            if self.old_key == "Antigen":
                self.methods["Antigen-Note"](key, value)
            else:
                note_instance = int(key.split("[")[1][:-1])
                note_index = [self.antibody_ann_dict[self.old_key].index(inst)
                              for inst in self.antibody_ann_dict[self.old_key]
                              if note_instance in inst["Instance"]][0]
                self.antibody_ann_dict[self.old_key][note_index]["Note"] = value.strip()
        else:
            self.antibody_ann_dict[self.old_key + "-" + key] = value.strip()
        # for item in self.antibody_ann_dict:
        #     print(f"Antibody Dict Item: {item}")

    def antigen_note(self, key, value):
        note_region = key[5:-1]
        for i in range(len(self.antibody_ann_dict[self.old_key])):
            if "Note" in self.antibody_ann_dict[self.old_key][i]:
                pass
            else:
                # print(f"Antigen dict 2: {self.antibody_ann_dict[self.old_key][i]}")
                if note_region == self.antibody_ann_dict[self.old_key][i]["Instance"]:
                    # print(f"Regions equal.")
                    self.antibody_ann_dict[self.old_key][i]["Note"] = value.strip()

    def any_instance_record(self, record):
        key, value = record.split(":")
        value = value.strip().split(" ") if " " in value else value.strip()

        name_key, instance = key.split("[")
        instance = instance[:-1].strip()
        if "," in instance:
            instance = instance.split(",")
        elif "-" in instance:
            instance = instance.split("-")

        # print(f"Outside all: {instance}")
        if all(item.isdigit() for item in instance):
            # print(f"In all instances: {instance}")
            instance = [int(num) for num in instance]

        if name_key not in self.antibody_ann_dict:
            data = [{"Instance": instance, name_key: value}]
            self.antibody_ann_dict[name_key] = data
        else:
            data = {"Instance": instance, name_key: value}
            self.antibody_ann_dict[name_key].append(data)

    def cdr_record(self, record):
        # key, value = record.split(":")
        # value = value.strip()
        # sequence = value.split(" ")[0]
        # residue = value.split(" ")[1][1:-1]
        # self.antibody_ann_dict[key] = sequence
        # self.antibody_ann_dict[key + "-Range"] = residue
        key, value = record.split(":")
        # value = value.strip().split(" ") if " " in value else value.strip()
        sequence = value.strip().split(" ")[0].strip()
        residue = value.strip().split(" ")[1][1:-1]

        if "]" not in key:
            self.antibody_ann_dict[key] = [{"Sequence": sequence, "Range": residue}]
            return None

        name_key, instance = key.split("[")
        instance = instance[:-1].strip()
        if "," in instance:
            instance = instance.split(",")
        elif "-" in instance:
            instance = instance.split("-")

        # print(f"Outside all: {instance}")
        if all(item.isdigit() for item in instance):
            # print(f"In all instances: {instance}")
            instance = [int(num) for num in instance]

        if name_key not in self.antibody_ann_dict:
            data = [{"Instance": instance, "Sequence": sequence, "Range": residue}]
            self.antibody_ann_dict[name_key] = data
        else:
            data = {"Instance": instance, "Sequence": sequence, "Range": residue}
            self.antibody_ann_dict[name_key].append(data)

    def heavy_chain_record(self, record):
        capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_sequence = ""

        value = record.split(":")[1].strip()
        for symbol in value:
            chain_sequence += symbol if symbol in capital_letters else ""

        key = record.split(":")[0]
        if "[" in key:
            heavy_chain_instances = key[key.index("[") + 1:-1]
            if "-" in heavy_chain_instances:
                heavy_chain_instances = [int(instance) for instance in heavy_chain_instances.split("-")]
            else:
                heavy_chain_instances = [int(instance) for instance in heavy_chain_instances.split(",")]
            self.antibody_ann_dict["Heavy Chain"] = [{"Instance": heavy_chain_instances,
                                                      "Sequence": chain_sequence}]
        else:
            self.antibody_ann_dict["Heavy Chain"] = chain_sequence

    def light_chain_record(self, record):
        capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_sequence = ""
        value = record.split(":")[1].strip()
        for symbol in value:
            chain_sequence += symbol if symbol in capital_letters else ""

        key = record.split(":")[0]
        if "[" in key:
            light_chain_instances = key[key.index("[") + 1:-1]
            if "-" in light_chain_instances:
                light_chain_instances = [int(instance) for instance in light_chain_instances.split("-")]
            else:
                light_chain_instances = [int(instance) for instance in light_chain_instances.split(",")]
            self.antibody_ann_dict["Light Chain"] = [{"Instance": light_chain_instances,
                                                      "Sequence": chain_sequence}]
        else:
            self.antibody_ann_dict["Light Chain"] = chain_sequence

    def chain_record(self, record):
        capital_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_sequence = ""
        value = record.split(":")[1].strip()
        for symbol in value:
            chain_sequence += symbol if symbol in capital_letters else ""

        key = record.split(":")[0]
        if "[" in key:
            chain_instances = key[key.index("[") + 1:-1]
            if "-" in chain_instances:
                chain_instances = [int(instance) for instance in chain_instances.split("-")]
            else:
                chain_instances = [int(instance) for instance in chain_instances.split(",")]
            self.antibody_ann_dict["Chain"] = [{"Instance": chain_instances,
                                                "Sequence": chain_sequence}]
        else:
            self.antibody_ann_dict["Chain"] = chain_sequence

    def normal_record(self, record):
        # print(f"Split Result: {record}")
        key, value = record.split(":", 1)
        if "[" in key:
            self.any_instance_record(record)
            # key = key.split("[")[0]
        else:
            value = value.strip()
            self.antibody_ann_dict[key] = value
