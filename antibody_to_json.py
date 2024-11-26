#  Populate JSON file from text data
import json


class AntibodyTxtJSON:

    def __init__(self):
        self.test = None

    @staticmethod
    def txt_to_dict(filename):
        antibody_dict = {}
        hanging_key = None
        chain = ""
        with open(filename, "r") as f:
            for line in f.readlines():
                if len(line) > 1:
                    if line.startswith("Note"):
                        note_line = line.strip().split(": ", 1)
                        antibody_dict[f"{split_line[0]}-Note"] = note_line[1]
                        continue

                    if line.startswith("Antigen"):
                        antigen_line = line.strip().split(": ", 1)
                        antibody_dict[f"{antigen_line[0]}"] = \
                            antigen_line[1][:antigen_line[1].index("(") - 1]
                        antibody_dict["Antigen-Gene"] = antigen_line[1].split()[-1][1:-2]
                        continue

                    if line.startswith("CDR") and not line.startswith("CDRSource"):
                        cdr_line = line.strip().split(": ", 1)
                        antibody_dict[f"{cdr_line[0]}"] = \
                            cdr_line[1][:cdr_line[1].index("(") - 1]
                        antibody_dict[f"{cdr_line[0]}-Range"] = cdr_line[1].split()[-1][1:-2]
                        continue

                    split_line = line.strip().split(": ", 1)

                    if len(split_line) == 2:
                        antibody_dict[split_line[0]] = split_line[1].replace(";", "")
                    elif len(split_line) == 1:
                        if ":" in split_line[0]:
                            hanging_key = split_line[0][:-1]
                        else:
                            if split_line[0] == "//":
                                antibody_dict[hanging_key] = chain
                                hanging_key = None
                                chain = ""
                            else:
                                chain += "".join(split_line[0].split(" ")[:-1]).strip()

        return antibody_dict


def main():
    filename = "11678.txt"
    antibody_json = AntibodyTxtJSON()
    antibody_dict_j = antibody_json.txt_to_dict(filename)
    with open(f"{filename.split('.')[0]}.json", "w") as jf:
        json.dump(antibody_dict_j, jf, indent=4)


if __name__ == "__main__":
    main()

