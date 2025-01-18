import antibody_ann_to_json


class AntibodyAnnRunner:

    def __init__(self, path):
        self.path = path
        self.aa_run = antibody_ann_to_json.AntibodyToJSON(path)

    def run(self):
        if len(self.aa_run.files) > 0:
            # Process the first file or any available one
            self.aa_run.single_file_transfer(self.path + "/" + self.aa_run.files[1])
        else:
            print("No files found in the directory:", self.path)

    # def run(self):
    #     self.aa_run.single_file_transfer(self.path + "/" + self.aa_run.files[25])


def main():
    # Run main Script
    runner = AntibodyAnnRunner("samples")
    runner.run()


if __name__ == "__main__":
    main()  # Comment 2 23.12.2024
