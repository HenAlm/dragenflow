import csv
import errno
import json
import os
from pathlib import Path
import shutil
from typing import List


def custom_sort(val: str) -> float:

    rank = 0.0
    if len(str(val)) > 1:
        if val[0] == "N":
            rank = float(val[-1]) - 0.5
        if val[0] == "T":
            rank = float(val[-1])

    return rank


def load_json(file: str = "config.json") -> dict:
    with open(file) as jf:
        configs = json.load(jf)
    return configs


def get_ref(excel: dict, template: dict) -> str:
    k_ = excel["RefGenome"]
    ref = template["ref_parameters"]["RefGenome"][k_]
    return ref["ref-dir"]


def set_fileprefix(excel: dict) -> str:
    prefix = "samplename"
    return prefix


def set_rgid(excel: dict) -> str:
    sample_sheet_path = excel["file_path"]
    flow_cell_id = get_flow_cell(sample_sheet_path)
    return flow_cell_id


def set_rgism(excel: dict) -> str:
    rgism = excel["SampleID"]
    return rgism


def create_fastq_dir(excel: list) -> list:
    for row in excel:
        path = Path(row["file_path"]).absolute()
        new_path = path.parent / row["SampleID"]
        new_path.mkdir(exist_ok=True)
        row["fastq_dir"] = new_path
    return excel


def fastq_file(excel: dict, read_n: int, copy_file: bool = True) -> str:

    sample_name = excel["Sample_Name"]
    sample_number = excel["Index"]
    lane_number = excel.get("Lane")
    if lane_number:
        file_name = (
            f"{sample_name}_S{sample_number}_L00{lane_number}_R{read_n}_001.fastq.gz"
        )

    else:
        file_name = f"{sample_name}_S{sample_number}_R{read_n}_001.fastq.gz"
    if copy_file:
        copy_fast_q(excel, file_name)
    return file_name


def copy_fast_q(excel: dict, fastq_f: str) -> None:

    sample_sheet_path = Path(excel["file_path"]).absolute().parent
    path_to_fastq = (
        sample_sheet_path / "demultiplex" / excel["Sample_Project"] / fastq_f
    )
    destination_of_fastq = Path(excel["fastq_dir"])
    if path_to_fastq.exists():
        shutil.copy(path_to_fastq, destination_of_fastq)
    else:
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), str(path_to_fastq)
        )


def check_key(dct: dict, k: str, val: str) -> dict:
    if k in dct.keys():
        dct[k] = val
        return dct
    else:
        raise KeyError


def dragen_cli(cmd: dict, excel: dict) -> str:
    default_str = " ".join(f"--{key} {val}" for (key, val) in cmd.items())
    final_str = f"grun.py -n dragen-{excel['Sample_Name']} -L logs -q  dragen.q -c  \
        dragen '{default_str}'"
    return final_str


def infer_pipeline(pipeline: str) -> str:
    str_list = pipeline.split("_")
    return str_list[0]


def get_flow_cell(path: str) -> str:
    split_path = path.split("/")
    flow_cell = split_path[-4]
    flow_cell_id = flow_cell.split("_")[-1]
    return flow_cell_id


def basic_reader(path: str) -> list:

    with open(path, newline="", encoding="utf-8") as inf:
        reader = csv.DictReader(inf)

        return list(reader)


def file_parse(path: str) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as inf:
        reader = csv.reader(inf)
        # find header row
        for row in reader:
            if "" not in row:
                fieldnames = row
                break
        else:
            # oops, *only* rows with empty cells found
            raise ValueError("Unable to determine header row")

        # rewind, switch to DictReader, skip past header
        # inf.seek(0)
        reader = csv.DictReader(inf, fieldnames)
        sorted_dict = sorted(
            reader, key=lambda row: custom_sort(row["tumor/normal"]), reverse=False
        )
        for dict_ in sorted_dict:
            dict_["file_path"] = path

        return sorted_dict