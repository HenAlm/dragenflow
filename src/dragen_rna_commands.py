import copy
from src.utility.dragen_utility import (
    fastq_file,
    set_fileprefix,
    set_rgid,
    set_rgism,
    get_ref_parameter,
)
from .utility.commands import Commands


class BaseDragenRnaCommand(Commands):
    """
    Base Dragen rna comands
    """

    def __init__(self, excel: dict, template: dict, seq_pipeline: str) -> None:
        self.excel = excel
        self.template = template
        self.seq_pipeline = seq_pipeline
        self.arg_registry = {
            "intermediate-results-dir": "/staging/intermediate",
            "output-file-prefix": set_fileprefix(self.excel),
            "tumor-fastq1": fastq_file(self.excel, 1),
            "tumor-fastq2": fastq_file(self.excel, 2),
            "RGID-tumor": set_rgid(self.excel),
            "RGSM-tumor": set_rgism(self.excel)
        }

    def construct_commands(self) -> dict:
        # get the arg from json config filie
        cmd_dict = copy.deepcopy(self.template[self.seq_pipeline])
        # get the dict that needs to be filled in at runtime
        param_list = [i for i in cmd_dict if str(cmd_dict[i]).startswith("{")]
        if len(param_list) == 0:
            raise RuntimeError("Someting went wrong with parsing template")
        for val in param_list:
            if val in self.arg_registry:
                cmd_dict[val] = self.arg_registry.get(val)
            else:
                tmp = cmd_dict[val][1:-1]
                cmd_dict[val] = get_ref_parameter(self.excel,self.template,tmp)
            if cmd_dict[val] == "":
                print(f"missing key '{val}' in registry or '{cmd_dict[val]}' in ref_parameters")
                continue
        return cmd_dict


class ExtraRnaCommands(Commands):
    """
    Custom Rna command specific certain instane
    """

    def construct_commands(self) -> dict:
        extra_cmd_dict = {}
        return extra_cmd_dict
