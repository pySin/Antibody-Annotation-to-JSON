#  Populate JSON file from text data


class AntibodyTxtJSON:

    def __init__(self):
        self.test = None

    @staticmethod
    def txt_to_dict(filename):
        hanging_key = None
        with open(filename, "r") as f:
            for line in f.readlines():
                if len(line) > 1:
                    split_line = line.strip().split(": ", 1)
                    if len(split_line) == 1:
                        hanging_key = split_line[0]
                    print(split_line)


def main():
    antibody_json = AntibodyTxtJSON()
    antibody_json.txt_to_dict("11678.txt")


if __name__ == "__main__":
    main()

