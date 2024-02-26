import json
import pathlib
import pandas as pd
from shutil import copyfile
import datetime

from qtutils import warning_dialog


def make_excel_backup(fn_excel: pathlib.Path) -> bool:

    ret = False
    if fn_excel.exists():
        fn_excel_bckup = (
            fn_excel.parent
            / "excel_backup"
            / (
                fn_excel.stem
                + "_"
                + datetime.datetime.now().strftime("%d%b%Y_%H%M%S")
                + ".xlsx"
            )
        )
        fn_excel_bckup.parents[0].mkdir(parents=True, exist_ok=True)
        ret = copyfile(fn_excel, fn_excel_bckup)
        if len(str(ret)) > 0:
            ret = True
    return ret


def return_nonempty_mnova_datasets(data: dict) -> dict:
    # remove datasets where multiplet counts is zero and peaks counts is zero and integrals counts is zero

    dicts_to_keep = {}
    for k, v in data.items():
        if isinstance(v, dict):
            if (
                (v["multiplets"]["count"] > 0)
                or (v["peaks"]["count"] > 0)
                or (v["integrals"]["count"] > 0)
            ):
                dicts_to_keep[k] = v
        else:
            dicts_to_keep[
                k
            ] = v  # must be smiles string, maybe need to check for this at some point

    return dicts_to_keep


def add_technique(key: str, technique_keys: dict, technique_counts: dict, expt: str):

    if expt in technique_keys.values():
        technique_counts[expt] += 1
        technique_keys[key] = f"{expt}_{technique_counts[expt]}"
    else:
        technique_keys[key] = expt


def read_in_mesrenova_json(fn: pathlib.Path) -> dict:
    """Read in the JSON file exported from MestReNova."""

    with open(fn, "r") as file:
        data_orig = json.load(file)

    data = return_nonempty_mnova_datasets(data_orig)

    # Identify the technique keys present in the JSON data
    technique_keys = {}
    technique_counts = {
        "HSQC": 0,
        "HMBC": 0,
        "COSY": 0,
        "C13_1D": 0,
        "H1_1D": 0,
        "NOESY": 0,
        "H1_pureshift": 0,
    }
    for key in data:
        if isinstance(data[key], dict):

            subtype = data[key].get("subtype", "")
            pulse_sequence = data[key].get("pulsesequence", "")

            if (
                subtype.lower().find("hsqc") != -1
                and data[key].get("type", "").lower() == "2d"
            ):
                add_technique(key, technique_keys, technique_counts, "HSQC")

            elif (
                subtype.lower().find("hmbc") != -1
                and data[key].get("type", "").lower() == "2d"
            ):
                add_technique(key, technique_keys, technique_counts, "HMBC")
            elif (
                subtype.lower().find("cosy") != -1
                and data[key].get("type", "").lower() == "2d"
            ):
                add_technique(key, technique_keys, technique_counts, "COSY")
            elif (
                subtype.lower().find("13c") != -1
                and data[key].get("type", "").lower() == "1d"
            ):
                add_technique(key, technique_keys, technique_counts, "C13_1D")
            elif (
                subtype.lower().find("1h") != -1
                and data[key].get("type", "").lower() == "1d"
                and "psyche" in data[key].get("pulsesequence", "").lower()
            ):
                add_technique(key, technique_keys, technique_counts, "H1_pureshift")
            elif (
                subtype.lower().find("1h") != -1
                and data[key].get("type", "").lower() == "1d"
            ):
                add_technique(key, technique_keys, technique_counts, "H1_1D")
            elif (
                subtype.lower().find("1h") != -1
                and data[key].get("type", "").lower() == "2d"
            ):
                add_technique(key, technique_keys, technique_counts, "NOESY")

    for k, v in technique_keys.items():
        data[v] = data[k]
        del data[k]

    return data


def get_2D_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    """
    Returns a pandas dataframe from the json_data dictionary for the specified technique.
    """
    df_data = []
    for i in range(json_data[technique]["peaks"]["count"]):
        df_data.append(
            [
                json_data[technique]["peaks"][str(i)]["delta2"],
                json_data[technique]["peaks"][str(i)]["delta1"],
                json_data[technique]["peaks"][str(i)]["intensity"],
                json_data[technique]["peaks"][str(i)]["type"],
            ]
        )

    df = pd.DataFrame(df_data, columns=["f2 (ppm)", "f1 (ppm)", "Intensity", "Type"])

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["f2 (ppm)"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df


def get_1d_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    df_data = []
    if json_data[technique]["multiplets"]["count"] == 0:
        # find peaks from  from peaks key
        for i in range(json_data[technique]["peaks"]["count"]):
            if str(i) in json_data[technique]["peaks"]:
                df_data.append(
                    [
                        json_data[technique]["peaks"][str(i)]["delta1"],
                        json_data[technique]["peaks"][str(i)]["intensity"],
                        json_data[technique]["peaks"][str(i)]["type"],
                    ]
                )

        df = pd.DataFrame(df_data, columns=["ppm", "Intensity", "Type"])

    else:
        # find peaks from  from multiplets key
        # Name	Shift	Range	H's	Integral	Class	J's	Method

        count = json_data[technique]["multiplets"]["count"]
        normValue = json_data[technique]["multiplets"]["normValue"]
        for i in [str(i) for i in range(count)]:
            if str(i) in json_data[technique]["multiplets"]:
                row = [
                    json_data[technique]["multiplets"][i]["delta1"],
                    json_data[technique]["multiplets"][i]["integralValue"],
                    json_data[technique]["multiplets"][i]["nH"],
                    json_data[technique]["multiplets"][i]["category"],
                ]

                # create a string from the list of J values and add it to df_data
                j_values = json_data[technique]["multiplets"][i]["jvals"]
                j_string = ", ".join([f"{j:1.3}" for j in j_values])
                j_string = f"{j_string}"
                row.append(j_string)
                df_data.append(row)

        df = pd.DataFrame(df_data, columns=["ppm", "Integral", "H's", "Class", "J's"])
        df["Integral"] = df["Integral"] / normValue

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df


def create_dataframes_from_mresnova_json(data: dict) -> dict:
    """
    Returns a dictionary of pandas dataframes for each technique in the data dictionary.
    """
    dataframes = {}
    for k, v in data.items():
        if k in [
            "H1_1D",
            "C13_1D",
            "HSQC",
            "HMBC",
            "HSQC_CH",
            "COSY",
            "NOESY",
            "H1_pureshift",
        ]:
            if v["type"].lower() == "2d":
                df = get_2D_dataframe_from_json(data, k)
                dataframes[k] = df
            elif v["type"].lower() == "1d":
                df = get_1d_dataframe_from_json(data, k)
                dataframes[k] = df
        elif k in ["smiles"]:
            dataframes["molecule"] = pd.DataFrame([data["smiles"]], columns=["smiles"])

    return dataframes


def write_excel_file_from_mresnova_df_dict(
    df_frames: dict, excel_path: pathlib.Path, qtstarted: bool = False
) -> bool:

    # check if path is valid
    if not excel_path.parent.exists():
        return False
    else:
        try:
            with pd.ExcelWriter(excel_path) as writer:
                # check if path is valid
                for k, df in df_frames.items():
                    df.to_excel(writer, sheet_name=k)
        except:
            # print exception
            warning_dialog(f"Exception occurred attempting to write excel file\n{str(excel_path)}",
                           "Exception occurred attempting to write excel file", qtstarted)
            return False
        return True


def check_for_multiple_HSQC_expts(data: dict) -> dict:
    """Check for multiple HSQC experiments and rename them."""

    hsqc_expts = [e for e in data.keys() if e.find("HSQC") != -1]
    print(hsqc_expts)

    if len(hsqc_expts) == 2:
        # compare the number of peaks in each experiment
        if (
            data[hsqc_expts[0]]["peaks"]["count"]
            >= data[hsqc_expts[1]]["peaks"]["count"]
        ):
            # rename the 2nd experiment
            data["HSQC_CH"] = data.pop(hsqc_expts[1])
        else:
            # rename the 1st experiment
            data["HSQC_CH"] = data.pop(hsqc_expts[0])
            # rename the 2nd experiment
            data["HSQC"] = data.pop(hsqc_expts[1])

    return data

    if __name__ == "__main__":
        pass
