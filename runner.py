import antibody_ann_to_json
import os


class AntibodyAnnRunner:

    def __init__(self, path):
        self.path = path
        self.aa_run = antibody_ann_to_json.AntibodyToJSON(path)

    @staticmethod
    def directory_json_files():
        directory = "json_files"
        os.makedirs(directory, exist_ok=True)

    def run(self):
        if len(self.aa_run.files) > 0:
            # Process the first file or any available one
            for i in range(len(self.aa_run.files)):
                print(f"Current File: {self.aa_run.files[i]}")
                print(f"File Number: {i}")
                self.aa_run.single_file_transfer(self.path + "/" + self.aa_run.files[i])
        else:
            print("No files found in the directory:", self.path)

    # def run(self):
    #     self.aa_run.single_file_transfer(self.path + "/" + self.aa_run.files[25])


def main():
    # Run main Script
    runner = AntibodyAnnRunner("samples")
    runner.directory_json_files()
    runner.run()


if __name__ == "__main__":
    main()  # Comment 2 23.12.2024
