#  Populate JSON file from text data


class AntibodyTxtJSON:

    def __init__(self):
        self.test = None

    @staticmethod
    def txt_to_dict(filename):
        hanging_key = None
        accumulate = []
        with open(filename, "r") as f:
            for line in f.readlines():
                if len(line) > 1:
                    split_line = line.strip().split(": ", 1)
                    if len(split_line) == 1 and ":" in split_line[0]:
                        hanging_key = split_line[0]
                        print(f"Inner key: {hanging_key}")
                    elif len(split_line) == 1:
                        accumulate.append(split_line[0])
                        print(f"Accumulate: {accumulate}")
                    else:
                        if hanging_key is not None:
                            print([hanging_key, accumulate])
                            hanging_key = None
                            accumulate = []
                        else:
                            print(split_line)
            if hanging_key is not None:
                print([hanging_key, accumulate])
                hanging_key = None
                accumulate = []


def main():
    antibody_json = AntibodyTxtJSON()
    antibody_json.txt_to_dict("11678.txt")


if __name__ == "__main__":
    main()

