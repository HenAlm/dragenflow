import logging
from typing import List, Optional

from .dragen_rna_commands import BaseDragenRnaCommand
from .utility.dragen_utility import (
    adapter_trimming,
    add_options,
    add_samplesheet_cols,
    dragen_cli,
    load_json,
    script_path, 
    SH_OVERRIDE,
)
from .utility.flow import Flow


class ConstructRnaPipeline(Flow):
    def constructor(self, excel: dict) -> Optional[List[str]]:
        self.profile = load_json(script_path("dragen_rna.json"))
        logging.info("executing dragen rna command")
        scripts = self.profile.get("scripts")
        if excel.get("disable_scripts"):
            scripts = None
        if self.profile["samplesheet"]:
            add_samplesheet_cols(excel,self.profile["samplesheet"])
        cmd_base = BaseDragenRnaCommand(excel, self.profile, "rna")
        cmd = cmd_base.construct_commands()
        cmd.update(adapter_trimming(self.profile, excel, cmd.get("read-trimmers")))
        cmd.update(add_options(excel[SH_OVERRIDE]))
        final_str = dragen_cli(cmd=cmd, excel=excel, scripts=scripts)
        return [final_str]
