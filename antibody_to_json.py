#  Populate JSON file from text data
import json


class AntibodyTxtJSON:

    def __init__(self):
        self.test = None

    @staticmethod
    def txt_to_dict(filename):
        antibody_dict = {}
        hanging_key = None
        accumulate = []
        with open(filename, "r") as f:
            for line in f.readlines():
                if len(line) > 1:
                    split_line = line.strip().split(": ", 1)

                    if len(split_line) == 2:
                        antibody_dict[split_line[0]] = split_line[1].replace(";", "")
                    elif len(split_line) == 1:
                        if ":" in split_line[0]:
                            hanging_key = split_line[0][:-1]
                            print(f"One value: {split_line}")
                        else:
                            if split_line[0] == "//":
                                antibody_dict[hanging_key] = accumulate
                                hanging_key = None
                                accumulate = []
                            else:
                                accumulate.append(split_line[0])

        #             if len(split_line) == 1 and ":" in split_line[0]:
        #                 hanging_key = split_line[0]
        #                 print(f"Inner key: {hanging_key}")
        #             elif len(split_line) == 1:
        #                 accumulate.append(split_line[0])
        #                 print(f"Accumulate: {accumulate}")
        #             else:
        #                 if hanging_key is not None:
        #                     antibody_dict[hanging_key] = accumulate
        #                     print([hanging_key, accumulate])
        #                     hanging_key = None
        #                     accumulate = []
        #                 else:
        #                     antibody_dict[split_line[0]] = split_line[1]
        #                     print(split_line)
        #         if hanging_key is not None:
        #             antibody_dict[hanging_key] = accumulate
        #             print([hanging_key, accumulate])
        #             hanging_key = None
        #             accumulate = []
        print(f"Antibody Dict: {antibody_dict}")
        return antibody_dict


def main():
    antibody_json = AntibodyTxtJSON()
    antibody_dict_j = antibody_json.txt_to_dict("11948.txt")
    with open("data.json", "w") as jf:
        json.dump(antibody_dict_j, jf, indent=4)


if __name__ == "__main__":
    main()

