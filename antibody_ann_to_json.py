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
            "Instance_Note": self.instance_note,
            "Multiple_Instance_Note": self.multiple_instance_note,
            "CDR": self.cdr_record,
            "Heavy Chain": self.heavy_chain_record,
            "Light Chain": self.light_chain_record,
            "Chain": self.chain_record,
            "Domains": self.domains_record,
            "Range": self.range_record,
            "MutationH": self.mutation_h_record,
            "MutationL": self.mutation_l_record,
            "HeavyPotentialNGlycos": self.heavy_potential_n_glycos_record,
            "HeavyConfirmedNGlycos": self.heavy_confirmed_n_glycos_record,
            "LightPotentialNGlycos": self.light_potential_n_glycos_record,
            "LightConfirmedNGlycos": self.light_confirmed_n_glycos_record,
            "LVGermline": self.lv_germline,
            "ConfirmedPTM": self.confirmed_ptm,
            "DisulfidesInter": self.disulfides_inter,
            "Linker": self.linker
        }

    @staticmethod
    def read_divide_records(filename):
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

            if "Format:" in records[0]:
                split_first_record = records[0].split("F")
                split_first_record[1] = "F" + split_first_record[1]
                records = records[1:]
                for i in range(len(split_first_record) -1, -1, -1):
                    records.insert(0, split_first_record[i])
            return records

    def single_file_transfer(self, filename):
        key_parts = ["Note", "CDR", "Range", "ConfirmedPTM", "DisulfidesInter"]  # Next:
        is_key_part_found = False

        # Produce JSON file from .txt annotation file
        current_records = self.read_divide_records(filename)
        for record in current_records:

            key = record.split(":")[0]
            key = key[:key.index("[")] if "[" in key else key

            for k_part in key_parts:
                if k_part in key:
                    self.methods[k_part](record)
                    is_key_part_found = True
                    break

            if is_key_part_found:
                is_key_part_found = False
                continue

            if key in self.methods:
                self.methods[key](record)
            else:
                self.normal_record(record)
            self.old_key = key

        filename = filename.split("/")[1]
        with open(f"{filename.split('.')[0]}.json", "w") as jf:
            json.dump(self.antibody_ann_dict, jf, indent=4)

    def antigen_record(self, record):
        key, value = record.split(":", 1)

        if "[" in key:
            name_key, instance = key.split("[")
            instance = [int(ins) for ins in instance[:-1].strip().split(",")]
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

    def range_record(self, record):
        key, value = record.split(":", 1)
        instance = int(key.split("[")[1][:-1])
        key = key.split("[")[0]
        start, end = value.strip().split(" ", 1)[0].split("-")

        mutations = value.split("(", 1)[1][:-2].split(" ") if "(" in value else "NONE"
        data = {"Instance": [instance], "Start": int(start), "End": int(end), "Mutations": mutations}

        if key in self.antibody_ann_dict:
            self.antibody_ann_dict[key].append(data)
        else:
            self.antibody_ann_dict[key] = [{"Instance": [instance], "Start": int(start), "End": int(end),
                                            "Mutations": mutations}]

    def domains_record(self, record):
        value = record.split(":", 1)[1].strip()

        if "Domains" not in self.antibody_ann_dict:
            self.antibody_ann_dict["Domains"] = ["dummy", value]
        else:
            self.antibody_ann_dict["Domains"].append(value)

    def note_record(self, record):
        key, value = record.split(":", 1)

        if key != "Note":
            if self.old_key == "Antigen":
                self.methods["Antigen-Note"](key, value)
            else:
                # note_instance = int(key.split("[")[1][:-1])
                # note_index = [self.antibody_ann_dict[self.old_key].index(inst)
                #               for inst in self.antibody_ann_dict[self.old_key]
                #               if note_instance in inst["Instance"]][0]
                # self.antibody_ann_dict[self.old_key][note_index]["Note"] = value.strip()
                # self.instance_note(key, value)
                self.methods["Instance_Note"](key, value)
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

    def instance_note(self, key, value):
        note_instance = key.split("[")[1][:-1]
        if "," in note_instance:
            self.methods["Multiple_Instance_Note"](key, value)
            return None

        note_instance = int(note_instance)
        note_index = [self.antibody_ann_dict[self.old_key].index(inst)
                      for inst in self.antibody_ann_dict[self.old_key]
                      if note_instance in inst["Instance"]][0]
        self.antibody_ann_dict[self.old_key][note_index]["Note"] = value.strip()

    def multiple_instance_note(self, key, value):
        note_instances = list(map(int, key.split("[")[1][:-1].split(",")))

        if self.old_key + "_Instances_Note" not in self.antibody_ann_dict:
            self.antibody_ann_dict[self.old_key + "_Instances_Note"] = [{"Instance": note_instances,
                                                                         "Note": value.strip()}]
        else:
            self.antibody_ann_dict[self.old_key + "_Instances_Note"].append({"Instance": note_instances,
                                                                             "Note": value.strip()})

    def any_instance_record(self, record):
        key, value = record.split(":")
        # value = value.strip().split(" ") if " " in value else value.strip()
        value = value.strip()  # [:-1]
        if any(v for v in value if v not in "1234567890 ;"):
            pass
        else:
            value = list(map(int, value.split(" ")))

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
            data = [{"Instance": instance, "Values": value}]
            self.antibody_ann_dict[name_key] = data
        else:
            data = {"Instance": instance, "Values": value}
            self.antibody_ann_dict[name_key].append(data)

    def cdr_record(self, record):
        key, value = record.split(":")
        sequence = value.strip().split(" ")[0].strip()

        if key == "CDRSource":
            self.antibody_ann_dict[key] = value.strip()
            return None

        if "CDRSource" in key:
            instances = [int(ins) for ins in key.split("[")[1][:-1].split(",")]
            key = key.split("[")[0]
            self.antibody_ann_dict[key] = [{"Sequence": instances, "Value": value.strip()}]
            return None
        start, end = list(map(int, value.strip().split(" ")[1][1:-1].split("-")))


        if "[" not in key:
            self.antibody_ann_dict[key] = [{"Sequence": sequence, "Start": start, "End": end}]
            return None

        name_key, instance = key.split("[")
        instance = instance[:-1].strip()
        if "," in instance:
            instance = instance.split(",")
        elif "-" in instance:
            instance = instance.split("-")

        if all(item.isdigit() for item in instance):
            instance = [int(num) for num in instance]

        if name_key not in self.antibody_ann_dict:
            data = [{"Instance": instance, "Sequence": sequence, "Start": start, "End": end}]
            self.antibody_ann_dict[name_key] = data
        else:
            data = {"Instance": instance, "Sequence": sequence, "Start": start, "End": end}
            self.antibody_ann_dict[name_key].append(data)

    def mutation_h_record(self, record):
        key, value = record.split(":", 1)
        instance = int(key.split("[")[1][:-1])
        mutations = value.split("(", 1)[0].strip().split(" ")
        reason = value.split("(", 1)[1][:-1]
        mutations_reasons = [{"Mutation": m, "Reason": reason} for m in mutations]

        if "MutationH" not in self.antibody_ann_dict:
            self.antibody_ann_dict["MutationH"] = [{"Instance": [instance],
                                                    "Mutations": mutations_reasons}]
        else:
            for i in range(len(self.antibody_ann_dict["MutationH"])):
                if self.antibody_ann_dict["MutationH"][i]["Instance"] == [instance]:
                    for current_reason in mutations_reasons:
                        self.antibody_ann_dict["MutationH"][i]["Mutations"].append(current_reason)
                    return None
            self.antibody_ann_dict["MutationH"].append({"Instance": [instance],
                                                        "Mutations": mutations_reasons})

    def mutation_l_record(self, record):
        key, value = record.split(":", 1)
        instance = int(key.split("[")[1][:-1])
        mutations = value.split("(", 1)[0].strip().split(" ")
        reason = value.split("(", 1)[1][:-2]
        mutations_reasons = [{"Mutation": m, "Reason": reason} for m in mutations]

        if "MutationL" not in self.antibody_ann_dict:
            self.antibody_ann_dict["MutationL"] = [{"Instance": [instance],
                                                    "Mutations": mutations_reasons}]
        else:
            for i in range(len(self.antibody_ann_dict["MutationL"])):
                if self.antibody_ann_dict["MutationL"][i]["Instance"] == [instance]:
                    for current_reason in mutations_reasons:
                        self.antibody_ann_dict["MutationL"][i]["Mutations"].append(current_reason)
                    return None
            self.antibody_ann_dict["MutationL"].append({"Instance": [instance],
                                                        "Mutations": mutations_reasons})

    def heavy_potential_n_glycos_record(self, record):
        print(f"Heavy Potential Record: {record}")
        # instance = int(record.split(":")[0].split("[", 1)[1][:-1])
        # confirmed_or_potential = "Confirmed" if "Confirmed" in record.split("[", 1)[0] else "Potential"
        value = record.split(":", 1)[1].strip()

        instance = record.split(":")[0].strip()
        if "[" not in instance:
            if "HeavyNGlycos" not in self.antibody_ann_dict:
                self.antibody_ann_dict["HeavyNGlycos"] = [{"Instance": "NONE", "Potential": [value]}]
                return None
            else:
                self.antibody_ann_dict["HeavyNGlycos"].append({"Instance": "NONE", "Potential": [value]})
                return None
        else:
            instance = instance.split("[", 1)[1][:-1]

        if "HeavyNGlycos" not in self.antibody_ann_dict:
            self.antibody_ann_dict["HeavyNGlycos"] = [{"Instance": [instance], "Potential": [value]}]
        else:
            for i in range(len(self.antibody_ann_dict["HeavyNGlycos"])):
                if self.antibody_ann_dict["HeavyNGlycos"][i]["Instance"][0] == instance:
                    if "Potential" in self.antibody_ann_dict["HeavyNGlycos"][i]:
                        return None
                    else:
                        self.antibody_ann_dict["HeavyNGlycos"][i]["Potential"] = [value]
                        return None

            self.antibody_ann_dict["HeavyNGlycos"].append({"Instance": [instance],
                                                           "Potential": [value]})

    def heavy_confirmed_n_glycos_record(self, record):
        # instance = int(record.split(":")[0].split("[", 1)[1][:-1])
        # confirmed_or_potential = "Confirmed" if "Confirmed" in record.split("[", 1)[0] else "Potential"
        value = record.split(":", 1)[1].strip()

        instance = record.split(":")[0].strip()
        if "[" not in instance:
            if "HeavyNGlycos" not in self.antibody_ann_dict:
                self.antibody_ann_dict["HeavyNGlycos"] = [{"Instance": "NONE", "Confirmed": [value]}]
                return None
            else:
                self.antibody_ann_dict["HeavyNGlycos"].append({"Instance": "NONE", "Confirmed": [value]})
                return None
        else:
            instance = instance.split("[", 1)[1][:-1]

        if "HeavyNGlycos" not in self.antibody_ann_dict:
            self.antibody_ann_dict["HeavyNGlycos"] = [{"Instance": [instance], "Confirmed": [value]}]
        else:
            for i in range(len(self.antibody_ann_dict["HeavyNGlycos"])):
                if self.antibody_ann_dict["HeavyNGlycos"][i]["Instance"][0] == instance:
                    if "Confirmed" in self.antibody_ann_dict["HeavyNGlycos"][i]:
                        return None
                    else:
                        self.antibody_ann_dict["HeavyNGlycos"][i]["Confirmed"] = [value]
                        return None

            self.antibody_ann_dict["HeavyNGlycos"].append({"Instance": [instance],
                                                           "Confirmed": [value]})

    def light_potential_n_glycos_record(self, record):
        instance = int(record.split(":")[0].split("[", 1)[1][:-1])
        # confirmed_or_potential = "Confirmed" if "Confirmed" in record.split("[", 1)[0] else "Potential"
        value = record.split(":", 1)[1].strip()

        if "LightNGlycos" not in self.antibody_ann_dict:
            self.antibody_ann_dict["LightNGlycos"] = [{"Instance": [instance], "Potential": [value]}]
        else:
            for i in range(len(self.antibody_ann_dict["LightNGlycos"])):
                if self.antibody_ann_dict["LightNGlycos"][i]["Instance"][0] == instance:
                    if "Potential" in self.antibody_ann_dict["LightNGlycos"][i]:
                        return None
                    else:
                        self.antibody_ann_dict["LightNGlycos"][i]["Potential"] = [value]
                        return None

            self.antibody_ann_dict["LightNGlycos"].append({"Instance": [instance],
                                                           "Potential": [value]})

    def light_confirmed_n_glycos_record(self, record):
        instance = int(record.split(":")[0].split("[", 1)[1][:-1])
        # confirmed_or_potential = "Confirmed" if "Confirmed" in record.split("[", 1)[0] else "Potential"
        value = record.split(":", 1)[1].strip()

        if "LightNGlycos" not in self.antibody_ann_dict:
            self.antibody_ann_dict["LightNGlycos"] = [{"Instance": [instance], "Confirmed": [value]}]
        else:
            for i in range(len(self.antibody_ann_dict["LightNGlycos"])):
                if self.antibody_ann_dict["LightNGlycos"][i]["Instance"][0] == instance:
                    if "Confirmed" in self.antibody_ann_dict["LightNGlycos"][i]:
                        return None
                    else:
                        self.antibody_ann_dict["LightNGlycos"][i]["Confirmed"] = [value]
                        return None

            self.antibody_ann_dict["LightNGlycos"].append({"Instance": [instance],
                                                           "Confirmed": [value]})

    def lv_germline(self, record):  # LVGermline[2]: Homo sapiens IGKV3-11*01;
        key, value = record.split(":", 1)
        instance = int(key.split("[", 1)[1][:-1])
        key = key.split("[")[0]
        species, gene_id = value.strip().rsplit(" ", 1)

        if key not in self.antibody_ann_dict:
            self.antibody_ann_dict[key] = [{"Instance": [instance], "Species": species, "GeneID": gene_id}]
        else:
            self.antibody_ann_dict[key].append({"Instance": [instance],
                                                "Species": species, "GeneID": gene_id})

    def confirmed_ptm(self, record):  # LightConfirmedPTM[2]: glycation 148 182 189 (rare);
        key, value = record.split(":", 1)
        instance = list(map(int, (key.split("[", 1)[1][:-1].split(","))))
        key = key.split("[")[0]
        # print(f'm_type:{value.strip().split(" ")[0]}')
        m_type = value.strip().split(" ")[0]
        frequency = value.strip().split(" ")[-1][1:-1] if "(" in value else ""
        position = [int(num) for num in value.strip().split(" ") if num.isnumeric()]

        if key not in self.antibody_ann_dict:
            self.antibody_ann_dict[key] = [{"Instance": instance, "Modifications":
                                           [{"Type": m_type, "Position": position, "Frequency": frequency}]}]
        else:
            for i in range(len(self.antibody_ann_dict[key])):
                if self.antibody_ann_dict[key][i]["Instance"] == instance:
                    self.antibody_ann_dict[key][i]["Modifications"].append({"Type": m_type,
                                                                            "Position": position,
                                                                            "Frequency": frequency})
                    return None
            self.antibody_ann_dict[key].append({"Instance": instance, "Modifications": [{"Type": m_type,
                                                "Position": position, "Frequency": frequency}]})

    def disulfides_inter(self, record):  # DisulfidesInterH1H2[2,5]: 222-230 225-233 350-353;
        key, value = record.split(":", 1)
        instance_a, instance_b = list(map(int, key.split("[")[1][:-1].split(",")))
        value = value.strip().split(" ")
        connections = [{"A": int(c.split("-")[0]), "B": int(c.split("-")[1])} for c in value]
        key = key.split("[")[0]

        self.antibody_ann_dict[key] = [{"InstanceA": instance_a, "InstanceB": instance_b, "Bonds": connections}]

    def linker(self, record):  # Linker[1,2]: 44-60;
        key, value = record.split(":", 1)
        l_from, l_to = [int(num) for num in key.split("[", 1)[1][:-1].split(",")]
        start, end = [int(num) for num in value.strip().split("-")]

        if "Linker" in self.antibody_ann_dict:
            self.antibody_ann_dict["Linker"].append({"From": l_from, "To": l_to, "Start": start, "End": end})
        else:
            self.antibody_ann_dict["Linker"] = [{"From": l_from, "To": l_to, "Start": start, "End": end}]

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
            self.antibody_ann_dict["HeavyChain"] = [{"Instance": heavy_chain_instances,
                                                      "Sequence": chain_sequence}]
        else:
            self.antibody_ann_dict["HeavyChain"] = chain_sequence

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
            self.antibody_ann_dict["LightChain"] = [{"Instance": light_chain_instances,
                                                      "Sequence": chain_sequence}]
        else:
            self.antibody_ann_dict["LightChain"] = chain_sequence

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
        key, value = record.split(":", 1)
        if "[" in key:
            self.any_instance_record(record)
        else:
            value = value.strip()
            self.antibody_ann_dict[key] = value
