import logging
from typing import List, Optional

from .dragen_met_commands import BaseDragenMetCommand
from .utility.dragen_utility import (
    adapter_trimming,
    add_options,
    add_samplesheet_cols,
    check_target,
    dragen_cli,
    load_json,
    script_path, 
    SH_OVERRIDE,
    SH_PARAM,
)
from .utility.flow import Flow


class ConstructMetPipeline(Flow):
    def constructor(self, excel: dict) -> Optional[List[str]]:
        self.profile = load_json(script_path("dragen_met.json"))
        logging.info("executing dragen methylation command")
        scripts = self.profile.get("scripts")
        if excel.get("disable_scripts"):
            scripts = None
        if self.profile["samplesheet"]:
            add_samplesheet_cols(excel,self.profile["samplesheet"])
        # make sure target gets set if given a named target
        check_target(
            excel, self.profile["ref_parameters"]["target"][excel["RefGenome"]]
        )
        cmd_base = BaseDragenMetCommand(excel, self.profile, excel[SH_PARAM])
        cmd = cmd_base.construct_commands()
        cmd.update(adapter_trimming(self.profile, excel, cmd.get("read-trimmers")))
        cmd.update(add_options(excel[SH_OVERRIDE]))
        final_str = dragen_cli(cmd=cmd, excel=excel, scripts=scripts)
        return [final_str]
