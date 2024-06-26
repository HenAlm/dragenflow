import argparse
import logging
from typing import List

from src.utility.flow import FlowConstructor
from src.dragen_pipeline import ConstructDragenPipeline
from src.dragen_met_pipeline import ConstructMetPipeline
from src.dragen_rna_pipeline import ConstructRnaPipeline
from src.utility.dragen_utility import (
    basic_reader,
    check_has_run,
    create_fastq_dir,
    file_parse,
    run_type,
    sort_list,
    SH_PARAM,
)

# register flows/pipeline
available_pipeline = {
    "dragen_dna": ConstructDragenPipeline(),
    "dragen_rna": ConstructRnaPipeline(),
    "dragen_met": ConstructMetPipeline(),
}

logging.basicConfig(filename="app.log", filemode="w", level=logging.DEBUG)
logging.info("started new logging session")


class HandleFlow(object):
    """
    Interface to flow(pipeline) objects through command line
    """

    def parse_file(self, path: str, flow: str) -> List[dict]:
        """Read excel file(sample sheet)

        if flow/pipeline is is dragen sort based on col tumor/normal
        else just read the csv file. Mehtod returns list of dictionary
        with column name as key and value as row.
        """
        if flow == "dragen":
            data_file = file_parse(path)
            return data_file
        else:
            data_file = basic_reader(path)
            return data_file

    def execute_bash(
        self,
        path: str,
        pipeline: str = "dragen",
        bash_cmd: str = "echo",
        dry_run: bool = False,
        disable_scripts: bool = False,
    ) -> list:
        """
        Construct bash command as string and execute if dry_run is False

        This creates appropriate flow object from argument supplied from cli
        & invoke construct_flow method of flow object
        """
        logging.info(f"dry run mode: {dry_run}")
        outputs = []
        command_list = []
        data_file = self.parse_file(path, pipeline)
        logging.info("creating fastq directory")
        data_file = create_fastq_dir(data_file, dry_run=dry_run)
        logging.info("assigning runtype")
        data_file1 = run_type(data_file)
        data_file = sort_list(data_file1)
        # chosen_pipeline = available_pipeline[pipeline]
        # flow_context = FlowConstructor(chosen_pipeline)
        for data in data_file:
            if data["pipeline"].lower() == "dragen":
                if data[SH_PARAM].startswith("rna"):
                    pipeline = "dragen_rna"
                    logging.info("Preparing dragen rna pipeline")
                elif data[SH_PARAM].startswith("methylation"):
                    pipeline = "dragen_met"
                    logging.info("Preparing dragen methylation pipeline")
                else:
                    pipeline = "dragen_dna"
                    logging.info("Preparing dragen dna pipeline")
                chosen_pipeline = available_pipeline[pipeline]
                flow_context = FlowConstructor(chosen_pipeline)
                # skip if pipeline is not dragen
                # attach script to data
                data["disable_scripts"] = disable_scripts
                logging.info("Creating dragen commands")
                constructed_str = flow_context.construct_flow(data=data)
                # contruct the commands first before checking as in case of paired sample
                # this would allow normal sample to have done previously and still be used
                if check_has_run(data):
                    logging.info(f"Skipping {data['fastq_dir']} as already executed.")
                    continue
                # collect all executable command in a list
                logging.info(f"Input dict:{data}")
                for c in constructed_str:
                    logging.info(f"command:{c}")
                    command_list.append([str(data["fastq_dir"]), c])
        if dry_run:
            for path, str_command in command_list:
                print("chdir " + path)
                outputs.append(str_command)
                print(str_command)
                print("===========")
        else:
            logging.info("Executing commands:")
            for path, str_command in command_list:
                output, arglist = FlowConstructor.execute_flow(
                    command=str_command, base_cmd=bash_cmd, wd_path=path
                )
                logging.info(f"Executed command: {arglist}")
                logging.info(f"Return code: {output.returncode}")
                outputs.append([(output.returncode, output.stdout)])
        return outputs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="dragenflow",
        description="Given a samplesheet, turn it into dragen commands.",
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        action="store",
        required=True,
        nargs=1,
        help="Required: path to the samplesheet file",
    )
    parser.add_argument(
        "-ds",
        "--disable_script",
        default=False,
        action="store_true",
        help="Optional: disable execution of pre/post script, defaults to True",
    )
    parser.add_argument(
        "-d",
        "--dryrun",
        default=False,
        action="store_true",
        help="Optional: enable dryrun, defaults to False",
    )
    parser.add_argument(
        "-c",
        "--cmd",
        default=None,
        help="Optional: if need to run arbitrary bash command, default None",
    )
    args = parser.parse_args()
    handle = HandleFlow()
    handle.execute_bash(
        path=args.path[0],
        bash_cmd=args.cmd,
        dry_run=args.dryrun,
        disable_scripts=args.disable_script,
    )
