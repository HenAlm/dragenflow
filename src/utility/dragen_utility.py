import csv
import errno
import json
import os
from pathlib import Path
import re
import shutil
from typing import List, Optional
import logging

# values for the samplesheet columns, SH_ for ones in file, SHA_ for added constructs
SHA_INDEX = 'row_index'
SHA_NPATH = "_normal_sample_path"
SHA_RTYPE = "_run_type"
SHA_TRG_NAME = "_target_name"
SH_NORMAL = "matching_normal_sample"
SH_OVERRIDE = "override"
SH_PARAM = "pipeline_parameters"
SH_SAMPLE = "SampleID"
SH_SM_PROJ = 'Sample_Project'
SH_TARGET = "TargetRegions"
SH_TUMOR = "Is_this_tumor"

# additional options selection
OPT_T_ALIGN = 'TUMOR_ALIGNMENT'
OPT_T_ANALYSIS = 'TUMOR_ANALYSIS'
OPTH_OPT = 'dragen option'
OPTH_VALUE = 'option value'
OPTH_SPEC = 'option specifier'

def custom_sort(val: str) -> float:
    rank = 0.0
    if len(str(val)) > 1:
        if val[0] == "N":
            rank = float(val[-1]) - 0.5
        if val[0] == "T":
            rank = float(val[-1])
    return rank


def script_path(filename: str) -> str:
    """
    A convenience function to get the absolute path to a file in this
    This allows the file to be launched from any
    directory.
    """
    filepath = os.path.join(os.path.dirname(__file__))
    config_path = os.path.join(filepath, "..")
    return os.path.join(config_path, filename)


def load_json(file: str = "config.json") -> dict:
    with open(file) as jf:
        configs = json.load(jf)
    return configs


def get_ref(excel: dict, template: dict) -> str:
    k_ = excel["RefGenome"]
    ref = template["ref_parameters"]["RefGenome"][k_]
    return ref["ref-dir"]


def set_fileprefix(excel: dict) -> str:
    sample_id = excel[SH_SAMPLE] if excel.get(SH_SAMPLE) else excel["Sample_ID"]
    return sample_id


def set_rgid(excel: dict) -> str:
    sample_sheet_path = excel["file_path"]
    flow_cell_id = get_flow_cell(sample_sheet_path)
    if excel.get("Lane"):
        flow_cell_id = f"{flow_cell_id}-{excel.get('Lane')}"
    return f"{flow_cell_id}-{excel[SHA_INDEX]}"


def set_rgism(excel: dict) -> str:
    rgism = excel[SH_SAMPLE] if excel.get(SH_SAMPLE) else excel["Sample_ID"]
    return rgism


def get_ref_parameter(excel: dict, template: dict, parameter:str) -> str:
    ref = excel["RefGenome"]
    if parameter in template["ref_parameters"]["RefGenome"][ref]:
        return template["ref_parameters"]["RefGenome"][ref][parameter]
    return ""


def create_fastq_dir(excel: list, dry_run: bool = False) -> List[dict]:
    for row in excel:
        if row["pipeline"].lower() != "dragen":
            continue
        path = Path(row["file_path"]).absolute()
        sample_id = row[SH_SAMPLE] if row.get(SH_SAMPLE) else row["Sample_ID"]
        new_path = path.parent / row[SH_SM_PROJ] / sample_id
        if not dry_run:
            new_path.mkdir(exist_ok=True)
        row["fastq_dir"] = new_path
        row["dry_run"] = dry_run
    return excel


def fastq_file(excel: dict, read_n: int, copy_file: bool = True) -> str:
    sample_name = excel["Sample_Name"]
    sample_number = excel[SHA_INDEX]
    lane_number = excel.get("Lane")
    if lane_number:
        file_name = (
            f"{sample_name}_S{sample_number}_L00{lane_number}_R{read_n}_001.fastq.gz"
        )
    else:
        file_name = f"{sample_name}_S{sample_number}_R{read_n}_001.fastq.gz"
    if copy_file:
        move_fast_q(excel, file_name)
    return file_name


def move_fast_q(excel: dict, fastq_f: str) -> None:
    sample_sheet_path = Path(excel["file_path"]).absolute().parent
    path_to_fastq = sample_sheet_path / excel[SH_SM_PROJ] / fastq_f
    destination_of_fastq = Path(excel["fastq_dir"])
    final_fastq_path = destination_of_fastq / fastq_f
    if excel["dry_run"]:
        if not (path_to_fastq.exists() or final_fastq_path.exists()):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), str(path_to_fastq)
            )
        return
    if path_to_fastq.exists():
        # in case if file already exist in destination
        if not final_fastq_path.exists():
            shutil.move(str(path_to_fastq), str(destination_of_fastq))
        log_path = destination_of_fastq / "logs"
        if not log_path.exists():
            os.mkdir(str(log_path))
    elif not final_fastq_path.exists():
        print("")
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), str(path_to_fastq)
        )


def check_key(dct: dict, k: str, val: str) -> dict:
    if k in dct.keys():
        dct[k] = val
        return dct
    else:
        raise KeyError


def dragen_cli(
    cmd: dict, excel: dict, postf: str = "", scripts: Optional[dict] = None
) -> str:
    default_str = " ".join(f"--{key} {val}" for (key, val) in cmd.items())
    grun_name = f"dragen-{excel['Sample_Name']}"
    if postf:
        grun_name = f"dragen-{excel['Sample_Name']}-{postf}"
    dragen_cmd = f"dragen {default_str}"
    if scripts:
        dragen_cmd = f"{scripts['pre']}\ndragen {default_str}\n{scripts['post']}"
    final_str = f"grun.py -n {grun_name} -L logs -q dragen.q -c '{dragen_cmd}'" # noqa: E501, B950
    return final_str


def infer_pipeline(pipeline: str) -> str:
    str_list = pipeline.split("_")
    return str_list[0]


def get_flow_cell(path: str) -> str:
    path = os.path.abspath(path)
    split_path = path.split("/")
    flow_cell = split_path[-3]
    flow_cell_id = flow_cell.split("_")[-1]
    return flow_cell_id


def basic_reader(path: str) -> list:
    with open(path, newline="", encoding="utf-8") as inf:
        reader = csv.DictReader(inf)
        return list(reader)


def trim_options(excel: dict, template: dict) -> str:
    # check if adaptertrim exist in sample sheet
    if not excel.get("AdapterTrim"):
        return ""
    elif excel["AdapterTrim"] == "truseq" or excel["AdapterTrim"] == "nextera":
        return template["adapters"][excel["AdapterTrim"]]
    elif excel["AdapterTrim"].startswith("/"):
        return excel["AdapterTrim"]
    else:
        raise ValueError(
            f"Cannot retrieve adapters for value: {excel['AdapterTrim']}"
        )


def adapter_trimming(template: dict, excel: dict, read_trimmer: str) -> dict:
    # take in read trimmers used and add adapter trimming if valid
    trim = trim_options(excel, template)
    cmd = {}
    if trim:
        if read_trimmer:
            cmd["read-trimmers"] = read_trimmer + ",adapter"
        else:
            cmd["read-trimmers"] = "adapter"
        cmd["trim-adapter-read1"] = trim
        cmd["trim-adapter-read2"] = trim
    return cmd


def file_parse(path: str, head_identifier="[Data]") -> List[dict]:
    # change the variable name
    with open(path, newline="", encoding="utf-8") as inf:
        reader = csv.reader(inf)
        # find header row
        for row in reader:
            if row[0].startswith(head_identifier):
                # if "" not in row:
                row = reader.__next__()
                fieldnames = row
                break
        else:
            # oops, *only* rows with empty cells found
            raise ValueError("Unable to determine header row")

        # if need to rewind, switch to DictReader, skip past header: inf.seek(0)
        reader = csv.DictReader(inf, fieldnames, restkey="__colmess")
        # Add the index before sorting
        row_index = 1
        reader = [val for val in reader]
        for row in reader:
            # convert Sample_ID into SampleID
            if row.get("Sample_ID"):
                row[SH_SAMPLE] = row.pop("Sample_ID")
            # if mistakes in samplesheet
            if row.get("__colmess"):
                raise ValueError(f"Sample {row[SH_SAMPLE]} has more columns than header")
            row[SHA_INDEX] = row_index
            row["file_path"] = path
            row_index += 1
            row[SHA_NPATH] = ""
            if SH_NORMAL in row and row[SH_NORMAL] and row[SH_NORMAL].startswith('/'):
                row[SHA_NPATH] = row[SH_NORMAL].rstrip('/')
                row[SH_NORMAL] = os.path.basename(row[SHA_NPATH])

        # Remove the pipeline that is not dragen
        # This must of be done before setting index, we need to preserve index
        reader = [row for row in reader if row["pipeline"] == "dragen"]

        return reader


def run_type(excel: List[dict]) -> List[dict]:
    for dt in excel:
        if dt[SH_PARAM] == "rna" or dt[SH_PARAM].startswith("methylation"):
            dt[SHA_RTYPE] = ""
            continue
        if (
            dt[SH_TUMOR] == "0"
            or dt[SH_TUMOR] == ""
            or dt[SH_TUMOR].lower() == "no"
        ):
            dt[SHA_RTYPE] = "germline"
        elif len(dt[SH_TUMOR]) >= 1 and dt[SH_NORMAL] == "":
            dt[SHA_RTYPE] = "somatic_single"
            if dt[SH_TUMOR] == "liquid":
                raise RuntimeError(
                    f"Liquid mode needs normal sample in {dt[SH_SAMPLE]} at {dt[SHA_INDEX]}"
                )
        elif len(dt[SH_TUMOR]) >= 1 and dt[SH_NORMAL] != "":
            sample_id = dt[SH_SAMPLE]
            normal_id = dt[SH_NORMAL]            
            sample_project = dt[SH_SM_PROJ]
            if check_sample(excel, normal_id, sample_project, dt[SHA_NPATH]):
                dt[SHA_RTYPE] = "somatic_paired"
            else:
                if dt[SHA_NPATH]:
                    normal_id = f"{dt[SHA_NPATH]}"
                raise RuntimeError(
                    f"Normal sample {normal_id} doesn't exist for {sample_id} at {dt[SHA_INDEX]}"
                )
        else:
            raise ValueError(f"invalid entry at index {dt[SHA_INDEX]}")
    return excel


def check_target(excel: dict, targets: dict) -> None:
    if not excel[SH_TARGET]:
        if excel[SH_PARAM] is "exome" or excel[SH_PARAM].startswith("umi"):
            raise ValueError(
                f"No target defined for {excel[SH_PARAM]} type in {excel[SH_SAMPLE]}."
            )
        return
    if excel[SH_TARGET].startswith("/"):
        return
    if not excel[SH_TARGET] in targets:
        raise ValueError(f"{SH_TARGET} value {excel[SH_TARGET]} not in json.")
    real_target = targets[excel[SH_TARGET]]
    excel[SHA_TRG_NAME] = excel[SH_TARGET]
    excel[SH_TARGET] = real_target
    return


def check_sample(excel: List[dict], sample_id: str, sample_project: str, sample_dir:str) -> bool:
    # if given path, check that it exists and that there is bam file
    if sample_dir:
        if not os.path.isdir(sample_dir):
            return False
        pref = os.path.basename(sample_dir)
        if os.path.isfile(os.path.join(sample_dir, pref + ".bam")):
            return True
        return False
    # some implicit assumption here that needs to be rechecked
    for dt in excel:
        if dt[SH_SAMPLE] == sample_id and dt[SH_SM_PROJ] == sample_project:
            return True
    return False


def sort_list(excel: List[dict], sorting_col: str = SHA_RTYPE) -> List[dict]:
    sorted_list = sorted(
        excel,
        key=lambda row: 1 if row[sorting_col] == "germline" else 0,
        reverse=True,
    )
    return sorted_list


def add_samplesheet_cols(excel:dict, add_cols:list=None, text_dir:str=None) -> None:
    # add given columns as json file into dir
    data = dict()
    if text_dir is None:
        text_dir = os.path.join(excel["fastq_dir"], 'logs')
    outfile = os.path.join(text_dir, "samplesheet_text.json")
    for i in add_cols:
        # add to json cols: value
        if i in excel:
            data[i] = excel[i]
    # if not sample dir at this point, then dryrun ... don't do anything'
    if not os.path.isdir(excel["fastq_dir"]):
        return
    # write dir/text.json
    if not os.path.isdir(text_dir):
        os.mkdir(text_dir)
    with open(outfile, 'w') as fs:
        json.dump(data,fs,sort_keys=True)


def check_has_run(excel:dict) -> bool:
    # check if sample command has executed: first find jobfiles
    jobfiles = []
    logs = f"{excel['fastq_dir']}/logs"
    if not os.path.isdir(logs):
        return False
    for fn in os.listdir(logs):
        if fn.endswith('.job'):
            jobfiles.append(fn)
    if len(jobfiles) == 0:
        return False
    # then check prefixes from commands
    prefix = []
    for fn in jobfiles:
        with open(f"{logs}/{fn}", 'r') as jobf:
            for line in jobf.readlines():
                if not line.startswith('dragen'):
                    continue
                m = re.search('--output-file-prefix (\S+)',line)
                if not m:
                    continue
                prefix.append(m.group(1))
    # finally check that all [prefix]-replay.json files exists
    for i in prefix:
        if not os.path.isfile(f"{excel['fastq_dir']}/{i}-replay.json"):
            return False
    return True


def add_options(opt_file:str, opt_target:str=None) -> dict:
    # read additional/update options from external file
    ret = dict()
    if not opt_file:
        return ret
    with open(opt_file,'r') as optfs:
        reader = csv.DictReader(optfs, dialect='excel-tab')
        for row in reader:
            if OPTH_OPT in row and row[OPTH_SPEC] and row[OPTH_SPEC] != opt_target:
                continue
            ret[row[OPTH_OPT]] = row[OPTH_VALUE]
    return ret
