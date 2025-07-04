import pathlib
import pandas as pd


class ConfigManager(object):
    @staticmethod
    def parse_region_size(size_str):
        return int(size_str.replace("kb", "")) * int(1e3)

    def __init__(self, config):
        self.config = config
        paths_file = self._format_patient_path(config["paths_file"])
        self.paths_df = pd.read_csv(paths_file, sep="\t")
        self.paths_df = self.paths_df.rename(columns={"cell_id": "sample"})

    @property
    def add_chr_prefix(self):
        return self.config.get("add_chr_prefix", True)

    # Input params
    @property
    def bin_sizes(self):
        return self.config.get("bin_sizes", ["1000kb"])

    @property
    def chromosomes(self):
        chroms = [str(x) for x in self.config["chromosomes"]]
        if "autosomes" in self.config["chromosomes"]:
            chroms.remove("autosomes")
            chroms = [str(x) for x in range(1, 23)] + chroms
        for c in chroms:
            assert c in [str(x) for x in range(1, 23)] + ["X"]
        if self.add_chr_prefix:
            chroms = ["chr{}".format(c) for c in chroms]
        return chroms

    @property
    def map_file(self):
        return pathlib.Path(self.config["map_file"])

    @property
    def min_bqual(self):
        return int(self.config.get("min_bqual", 30))

    @property
    def min_mqual(self):
        return int(self.config.get("min_mqual", 30))

    @property
    def normal_bam_file(self):
        path = pathlib.Path(self.config["normal_bam_file"])
        path = self._format_patient_path(path)
        return path.absolute()

    @property
    def out_dir(self):
        path = pathlib.Path(self.config["out_dir"])
        path = self._format_patient_path(path)
        return path.absolute()

    @property
    def patient_id(self):
        return self.config.get("patient_id", None)

    @property
    def pipeline_dir(self):
        path = pathlib.Path(self.config["pipeline_dir"])
        path = self._format_patient_path(path)
        return path.absolute()

    @property
    def log_dir(self):
        return self.pipeline_dir.joinpath("log")

    @property
    def ref_dir(self):
        return self.working_dir.joinpath("ref")

    @property
    def ref_genome_file(self):
        path = pathlib.Path(self.config["ref_genome_file"])
        return path.absolute()

    @property
    def working_dir(self):
        return self.pipeline_dir.joinpath("working")

    # Templates
    @property
    def allele_bin_counts_template(self):
        return self.out_dir.joinpath("allele_bin_counts", "{prog}", "{bin_size}", "{sample}.tsv.gz")

    @property
    def allele_counts_chrom_template(self):
        return self.working_dir.joinpath("allele_counts_chrom", "{prog}", "{chrom}", "{sample}.tsv.gz")

    @property
    def allele_counts_template(self):
        return self.out_dir.joinpath("allele_counts", "{prog}", "{sample}.tsv.gz")

    @property
    def gc_template(self):
        return self.ref_dir.joinpath("gc_{bin_size}.tsv")

    @property
    def map_template(self):
        return self.ref_dir.joinpath("map_{bin_size}.tsv")

    @property
    def rdr_corrected_template(self):
        return self.out_dir.joinpath("rdr_corrected", "{bin_size}", "{sample}.tsv.gz")

    @property
    def rdr_template(self):
        return self.working_dir.joinpath("rdr", "{bin_size}", "{sample}.tsv.gz")

    @property
    def reads_chrom_template(self):
        return self.working_dir.joinpath("reads_chrom", "{bin_size}", "{chrom}", "{sample}.tsv.gz")

    @property
    def reads_template(self):
        return self.working_dir.joinpath("reads", "{bin_size}", "{sample}.tsv.gz")

    @property
    def baf_plot_template(self):
        return self.out_dir.joinpath("plots", "baf", "{prog}", "{bin_size}", "{sample}.png")

    @property
    def rdr_plot_template(self):
        return self.out_dir.joinpath("plots", "rdr", "{bin_size}", "{sample}.png")

    # File getters
    def get_pipeline_files(self):
        for s in self.paths_df["sample"].unique():
            for b in self.bin_sizes:
                yield str(self.rdr_corrected_template).format(bin_size=b, sample=s)

                yield str(self.rdr_plot_template).format(bin_size=b, sample=s)

        for s in self.paths_df["sample"].unique():
            for p in self.config["snp_files"]:
                for h in self.bin_sizes:
                    yield str(self.allele_bin_counts_template).format(bin_size=h, prog=p, sample=s)

                    yield str(self.baf_plot_template).format(bin_size=h, prog=p, sample=s)

    def get_bam_file(self, wc):
        if wc.sample == "normal":
            return self.normal_bam_file
        df = self.paths_df.set_index("sample")
        return df.loc[wc.sample]["path"]

    def get_snp_file(self, wc):
        return self._format_patient_path(self.config["snp_files"][wc.prog])

    def get_log_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)

    def _format_patient_path(self, path):
        if self.patient_id is not None:
            path = pathlib.Path(str(path).format(patient_id=self.patient_id))
        return path

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.working_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
