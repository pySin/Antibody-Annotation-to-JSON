import antibody_ann_to_json


class AntibodyAnnRunner:

    def __init__(self, path):
        self.path = path
        self.aa_run = antibody_ann_to_json.AntibodyToJSON(path)

    def run(self):
        self.aa_run.single_file_transfer(self.path + "/" + self.aa_run.files[3])


def main():
    runner = AntibodyAnnRunner("samples")
    runner.run()


if __name__ == "__main__":
    main()
