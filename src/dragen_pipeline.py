import logging
from typing import List, Optional

from .dragen_commands import (
    BaseDragenCommand,
    PairedVariantCommands,
)
from .utility.commands import CompositeCommands
from .utility.dragen_utility import (
    check_target,
    dragen_cli,
    load_json,
    script_path,
    SH_TARGET,
    SHA_RTYPE,
    SHA_TRG_NAME,
    trim_options,
)
from .utility.flow import Flow


class ConstructDragenPipeline(Flow):
    def __init__(self):
        self.normals = {}
        self.profile = None

    def check_trimming(self, excel: dict, read_trimmer) -> dict:
        trim = trim_options(excel, self.profile)
        cmd = {}
        cmd["read-trimmers"] = read_trimmer
        if trim:
            cmd["read-trimmers"] = cmd["read-trimmers"] + ",adapter"
            cmd["trim-adapter-read1"] = trim
            cmd["trim-adapter-read2"] = trim
        return cmd

    def add_cnv(self, excel: dict, cmd: dict) -> bool:
        tmp = self.profile["ref_parameters"]["cnvpanelofnormals"]
        if SHA_TRG_NAME not in excel or not excel[SHA_TRG_NAME]:
            return False
        if not excel[SHA_TRG_NAME] in tmp[excel["RefGenome"]]:
            return False
        cmd["cnv-normals-list"] = tmp[excel["RefGenome"]][excel[SHA_TRG_NAME]]
        cmd["cnv-target-bed"] = excel[SH_TARGET]
        cmd["enable-cnv"] = "true"
        return True

    def command_with_trim(self, excel: dict, pipe_elem: str) -> dict:
        pipeline = excel.get("pipeline_parameters")
        base_cmd = BaseDragenCommand(excel, self.profile, f"{pipeline}_{pipe_elem}")
        cmd = base_cmd.construct_commands()
        trim_cmd = self.check_trimming(excel, cmd.get("read-trimmers"))
        return {**cmd, **trim_cmd}

    def umi_pipeline(self, excel:dict, pipe_elem:str, tumor:bool=False) -> dict:
        pipeline = excel.get("pipeline_parameters")
        cmd_base = BaseDragenCommand(
            excel, self.profile, f"{pipeline}_{pipe_elem}"
        )
        cmd_base.set_umi_fastq(excel, tumor)
        return cmd_base.construct_commands()

    def sample_pon(self, key: str, dryrun: bool, sample_dir: str, cmd: dict) -> None:
        # create temporary cnv pon with normal added
        add_normal = f"{self.normals[key]}.target.counts.gc-corrected.gz"
        new_panel = f"{sample_dir}/logs/cnv_pon.txt"
        if not dryrun:
            with open(new_panel, 'w') as new_list:
                with open(cmd["cnv-normals-list"], 'r') as old_list:
                    for line in old_list.readlines():
                        new_list.write(line)
                new_list.write(add_normal)
        cmd["cnv-normals-list"] = new_panel

    def constructor(self, excel: dict) -> Optional[List[str]]:
        self.profile = load_json(script_path("dragen_config.json"))["profile1"]
        # load pre and post scripts
        scripts = self.profile.get("scripts")
        if excel.get("disable_scripts"):
            scripts = None
        # make sure target gets set if given a named target
        check_target(
            excel, self.profile["ref_parameters"]["target"][excel["RefGenome"]]
        )

        # no pipeline set, check if target to choose between exome and genome
        if not excel["pipeline_parameters"]:
            if excel[SH_TARGET]:
                excel["pipeline_parameters"] = "exome"
            else:
                excel["pipeline_parameters"] = "genome"
        pipeline = excel["pipeline_parameters"]

        if excel[SHA_RTYPE] == "germline":
            # out put prefix = samplename
            # out put prefix paired = sample.s
            if pipeline.startswith("umi"):
                logging.info(f"{excel[SHA_RTYPE]}: executing umi normal_pipeline")
                cmd_d = self.umi_pipeline(excel, "normal_pipeline")
            else:
                logging.info(f"{excel[SHA_RTYPE]}: executing normal_pipeline")
                cmd_d = self.command_with_trim(excel, "normal_pipeline")
                self.add_cnv(excel, cmd_d)
            # store bam file
            self.normals[
                f"{excel['Sample_Project']}/{excel['SampleID']}"
            ] = f"../{excel['SampleID']}/{cmd_d['output-file-prefix']}"
            final_str = dragen_cli(cmd=cmd_d, excel=excel, scripts=scripts)
            return [final_str]

        elif excel[SHA_RTYPE] == "somatic_single":
            if pipeline.startswith("umi"):
                logging.info(f"{excel[SHA_RTYPE]}: executing umi tumor_pipeline")
                cmd_d = self.umi_pipeline(excel, "tumor_pipeline", True)
            else:
                logging.info(f"{excel[SHA_RTYPE]}: executing tumor_pipeline")
                cmd_d = self.command_with_trim(excel, "tumor_pipeline")
                self.add_cnv(excel, cmd_d)
            final_str = dragen_cli(cmd=cmd_d, excel=excel, scripts=scripts)
            return [final_str]

        elif excel[SHA_RTYPE] == "somatic_paired":
            arg_string = []
            base_cmd = BaseDragenCommand(
                excel, self.profile, f"{pipeline}_paired_variant_call"
            )
            cmd_d2 = CompositeCommands()
            normal_prefix = (
                f"{excel['Sample_Project']}/{excel['matching_normal_sample']}"
            )
            if pipeline.startswith("umi"):
                # step 1
                logging.info(f"{excel[SHA_RTYPE]}: preparing umi alignment template")
                cmd_d1 = self.umi_pipeline(excel, "tumor_alignment", True)
                final_str1 = dragen_cli(
                    cmd=cmd_d1, excel=excel, postf="alignment", scripts=scripts
                )
                arg_string.append(final_str1)
            else:
                logging.info(f"{excel[SHA_RTYPE]}: preparing tumor alignment template")
                # step 1 tumor alignment
                cmd_d1 = self.command_with_trim(excel, "tumor_alignment")
                if self.add_cnv(excel, cmd_d1):
                    cmd_d1["enable-map-align-output"] = "true"
                    self.sample_pon(normal_prefix, excel["dry_run"], excel["fastq_dir"], cmd_d1)
                final_str1 = dragen_cli(
                    cmd=cmd_d1, excel=excel, postf="alignment", scripts=scripts
                )
                arg_string.append(final_str1)
            # step 2 paired variant call
            logging.info(f"{excel[SHA_RTYPE]}: preparing paired variant call template")
            pv_cmd = PairedVariantCommands(f"{self.normals[normal_prefix]}.bam", cmd_d1)
            pv_cmd.add_error_cal(excel[SH_TARGET])
            cmd_d2.add(base_cmd)
            cmd_d2.add(pv_cmd)
            final_str2 = dragen_cli(
                cmd=cmd_d2.construct_commands(),
                excel=excel,
                postf="analysis",
                scripts=scripts,
            )
            arg_string.append(final_str2)
            return arg_string

        else:
            logging.info(
                f"No known pipeline run type, problem on \
                    {excel['Sample_Project']}:{excel.get('SampleID')}, skipping"
            )
