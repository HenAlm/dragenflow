import copy

from .utility.commands import Commands
from .utility.dragen_utility import (
    fastq_file,
    set_fileprefix,
    set_rgid,
    set_rgism,
    get_ref_parameter,
    SH_PARAM,
    SH_TARGET,
)


class BaseDragenCommand(Commands):
    """
    Base Dragen commands
    """

    def __init__(self, excel: dict, template: dict, seq_pipeline: str) -> None:
        self.excel = excel
        self.template = template
        self.seq_pipeline = seq_pipeline
        self.arg_registry = {
            "fastq-file1": fastq_file(self.excel, 1),
            "fastq-file2": fastq_file(self.excel, 2),
            "tumor-fastq1": fastq_file(self.excel, 1),
            "tumor-fastq2": fastq_file(self.excel, 2),
            "output-file-prefix": set_fileprefix(self.excel),
            "qc-coverage-region-1": self.excel[SH_TARGET],
            "RGID": set_rgid(self.excel),
            "RGSM": set_rgism(self.excel),
            "RGID-tumor": set_rgid(self.excel),
            "RGSM-tumor": set_rgism(self.excel),
            # depending on the use case this can be directly added to json-template file
            "intermediate-results-dir": "/staging/intermediate",
            "vc-snp-error-cal-bed": self.excel[SH_TARGET]
        }

    def construct_commands(self) -> dict:
        # select the parameter from config template
        cmd_dict1 = self.template[self.seq_pipeline]
        cmd_dict2 = copy.deepcopy(cmd_dict1)
        # get the dict that needs to be filled in at runtime
        param_list = [i for i in cmd_dict2 if str(cmd_dict2[i]).startswith("{")]
        if len(param_list) == 0:
            raise RuntimeError("Something wrong with parsing template")
        for val in param_list:
            if val in self.arg_registry:
                cmd_dict2[val] = self.arg_registry.get(val)
            else:
                tmp = cmd_dict2[val][1:-1]
                cmd_dict2[val] = get_ref_parameter(self.excel,self.template,tmp)
            if cmd_dict2[val] == "":
                print(f"missing key '{val}' in registry or '{cmd_dict2[val]}' in ref_parameters")
                continue
        return cmd_dict2

    def set_umi_fastq(self, excel: dict, is_tumor: bool = False) -> None:
        # if normal umis, need to swap fastqs around
        if excel[SH_PARAM] == "umi":
            if is_tumor:
                self.arg_registry["tumor-fastq2"] = fastq_file(self.excel, 3)
            else:
                self.arg_registry["fastq-file2"] = fastq_file(self.excel, 3)
            self.arg_registry["umi-fastq"] = fastq_file(self.excel, 2, False)
        self.arg_registry["umi-metrics-interval-file"] = excel[SH_TARGET]
        return


class PairedVariantCommands(Commands):
    """
    PairedVariant specific dragen command
    """

    def __init__(self, normal: str, tumor_align: dict, tumor_varc: dict) -> None:
        self.tumor_a = tumor_align
        self.tumor_vc = tumor_varc
        self.normal_pref = normal

    def construct_commands(self) -> dict:
        cmd_dict = {}
        cmd_dict["bam-input"] = self.tumor_vc["bam-input"].replace("{normal}",self.normal_pref)
        cmd_dict["tumor-bam-input"] = f"{self.tumor_a['output-file-prefix']}_tumor.bam"
        cmd_dict["output-file-prefix"] = f"{self.tumor_a['output-file-prefix']}.tn"
        if "cnv-normal-cnv-vcf" in self.tumor_vc:
            opts = ["cnv-normal-cnv-vcf","cnv-normal-b-allele-vcf"]
            for i in opts:
                cmd_dict[i] = self.tumor_vc[i].replace("{normal}",self.normal_pref)
        if "vc-snp-error-cal-bed" in self.tumor_vc:
            cmd_dict["vc-snp-error-cal-bed"] = self.tumor_a["qc-coverage-region-1"]
        return cmd_dict
