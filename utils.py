import pathlib


class MyConfigManager(object):
    @staticmethod
    def parse_region_size(size_str):
        return int(size_str.replace("kb", "")) * int(1e3)

    def __init__(self, config):
        self.config = config

    # Parameters
    @property
    def bin_size_baf(self):
        return self.config.get("bin_size_baf", "10kb")

    @property
    def bin_size_rdr(self):
        return self.config.get("bin_size_rdr", "100kb")

    @property
    def chromosomes(self):
        return [str(x) for x in self.config["chromosomes"]]

    @property
    def chromosome_parallel(self):
        return self.config.get("chromosome_parallel", True)

    @property
    def map_file(self):
        return pathlib.Path(self.config["map_file"])

    @property
    def min_bqual(self):
        return int(self.config.get("min_bqual", 30))

    @property
    def min_mqual(self):
        return int(self.config.get("min_mqual", 30))

    # Input files
    @property
    def normal_bam_file(self):
        return pathlib.Path(self.config["normal_bam_file"])

    @property
    def tumour_bam_files(self):
        return {k: pathlib.Path(v) for k, v in self.config["tumour_bam_files"].items()}

    @property
    def ref_genome_file(self):
        return pathlib.Path(self.config["ref_genome_file"])

    @property
    def snp_file(self):
        return pathlib.Path(self.config["snp_file"])

    # Directories
    @property
    def log_dir(self):
        return self.pipeline_dir.joinpath("log")

    @property
    def out_dir(self):
        return pathlib.Path(self.config["out_dir"])

    @property
    def pipeline_dir(self):
        return pathlib.Path(self.config["pipeline_dir"])

    @property
    def ref_dir(self):
        return self.working_dir.joinpath("ref")

    @property
    def working_dir(self):
        return self.pipeline_dir.joinpath("working")

    # Templates
    @property
    def allele_bin_counts_template(self):
        return self.out_dir.joinpath("allele_bin_counts", "{sample}.tsv.gz")

    @property
    def allele_counts_chrom_template(self):
        return self.working_dir.joinpath("allele_counts_chrom", "{chrom}", "{sample}.tsv.gz")

    @property
    def allele_counts_template(self):
        return self.out_dir.joinpath("allele_counts", "{sample}.tsv.gz")

    @property
    def gc_template(self):
        return self.ref_dir.joinpath("gc.tsv")

    @property
    def map_template(self):
        return self.ref_dir.joinpath("map.tsv")

    @property
    def rdr_corrected_template(self):
        return self.out_dir.joinpath("rdr_corrected", "{sample}.tsv.gz")

    @property
    def rdr_template(self):
        return self.working_dir.joinpath("rdr", "{sample}.tsv.gz")

    @property
    def reads_chrom_template(self):
        return self.working_dir.joinpath("reads_chrom", "{chrom}", "{sample}.tsv.gz")

    @property
    def reads_template(self):
        return self.working_dir.joinpath("reads", "{sample}.tsv.gz")

    @property
    def baf_plot_template(self):
        return self.out_dir.joinpath("plots", "baf", "{sample}.png")

    @property
    def rdr_plot_template(self):
        return self.out_dir.joinpath("plots", "rdr", "{sample}.png")

    # File getters
    def get_pipeline_files(self):
        for s in self.tumour_bam_files:
            yield str(self.rdr_corrected_template).format(sample=s)

            yield str(self.rdr_plot_template).format(sample=s)

            yield str(self.allele_bin_counts_template).format(sample=s)

            yield str(self.baf_plot_template).format(sample=s)

    def get_bam_file(self, wc):
        if wc.sample == "normal":
            return self.normal_bam_file
        else:
            return self.tumour_bam_files[wc.sample]

    def get_log_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.working_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
