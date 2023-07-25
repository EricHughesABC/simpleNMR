
import json
import pathlib
import pandas as pd


def return_nonempty_mnova_datasets(data: dict)->dict:
    # remove datasets where multiplet counts is zero and peaks counts is zero and integrals counts is zero

    dicts_to_keep = {}
    for k, v in data.items():
        if isinstance(v, dict):
            if( v['multiplets']["count"] > 0) or (v['peaks']["count"] > 0) or (v['integrals']["count"] > 0):
                dicts_to_keep[k] = v
        else:
            dicts_to_keep[k] = v  # must be smiles string, maybe need to check for this at some point

    return dicts_to_keep


def read_in_mesrenova_json(fn: pathlib.Path)->dict:
    """Read in the JSON file exported from MestReNova."""

    with open(fn, 'r') as file:
        data_orig = json.load(file)

    data = return_nonempty_mnova_datasets(data_orig)
    print("kept datasets = ", data.keys())

    # Identify the technique keys present in the JSON data
    technique_keys = {}
    for key in data:
        if isinstance(data[key], dict):

            
            subtype = data[key].get('subtype', '')
            pulse_sequence = data[key].get('pulsesequence', '')

            if subtype.lower().find('hsqc') != -1 and data[key].get('type', '').lower() == '2d':
                technique_keys[key] = 'HSQC'
            elif subtype.lower().find('hmbc') != -1 and data[key].get('type', '').lower() == '2d':
                technique_keys[key] = 'HMBC'
            elif subtype.lower().find('cosy') != -1 and data[key].get('type', '').lower() == '2d':
                technique_keys[key] = 'COSY'
            elif subtype.lower().find('13c') != -1 and data[key].get('type', '').lower() == '1d':
                technique_keys[key] = 'C13_1D'
            elif subtype.lower().find('1h') != -1 and data[key].get('type', '').lower() == '1d' and "psyche" in data[key].get('pulsesequence', '').lower():
                technique_keys[key] ='H1_pureshift'
            elif subtype.lower().find('1h') != -1 and data[key].get('type', '').lower() == '1d':
                technique_keys[key] ='H1_1D'
            elif subtype.lower().find('1h') != -1 and data[key].get('type', '').lower() == '2d':
                technique_keys[key] ='NOESY'

            if  pulse_sequence.lower().find('hmbc') != -1 and data[key].get('type', '').lower() == '2d':
                technique_keys[key] = 'HMBC'

    for k, v in technique_keys.items():
        data[v] = data[k]
        del data[k]

    return data

def get_2D_dataframe_from_json(json_data: dict, technique: str)->pd.DataFrame:
    """
    Returns a pandas dataframe from the json_data dictionary for the specified technique.
    """
    df_data = []
    for i in range(json_data[technique]["peaks"]["count"]):
        df_data.append([json_data[technique]["peaks"][str(i)]["delta2"], 
                        json_data[technique]["peaks"][str(i)]["delta1"], 
                        json_data[technique]["peaks"][str(i)]["intensity"], 
                        json_data[technique]["peaks"][str(i)]["type"]])

    df = pd.DataFrame(df_data, columns=["f2 (ppm)", "f1 (ppm)", "Intensity", "Type"])

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["f2 (ppm)"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df

def get_1d_dataframe_from_json( json_data: dict, technique: str)->pd.DataFrame:
    df_data = []
    if json_data[technique]["multiplets"]["count"] == 0:
        # find peaks from  from peaks key
        for i in range(json_data[technique]["peaks"]["count"]):
            if str(i) in json_data[technique]["peaks"]:
                df_data.append([json_data[technique]["peaks"][str(i)]["delta1"], 
                                json_data[technique]["peaks"][str(i)]["intensity"],
                                json_data[technique]["peaks"][str(i)]["type"]])
                
        df = pd.DataFrame(df_data, columns=["ppm",  "Intensity", "Type"])

    else:
        # find peaks from  from multiplets key
        # Name	Shift	Range	H's	Integral	Class	J's	Method

        count = json_data[technique]["multiplets"]["count"]
        normValue = json_data[technique]["multiplets"]["normValue"]
        for i in [str(i) for i in range(count)]:
            if str(i) in json_data[technique]["multiplets"]:
                row = [json_data[technique]["multiplets"][i]["delta1"], 
                                json_data[technique]["multiplets"][i]["integralValue"],
                                json_data[technique]["multiplets"][i]["nH"],
                                json_data[technique]["multiplets"][i]["category"]]
                
                # create a string from the list of J values and add it to df_data
                j_values = json_data[technique]["multiplets"][i]["jvals"]
                j_string = ", ".join([f"{j:1.3}" for j in j_values])
                j_string = f"{j_string}"
                row.append(j_string)
                df_data.append(row)

        
        df = pd.DataFrame(df_data, columns=["ppm", "Integral",  "H's", "Class", "J's"])
        df["Integral"] = df["Integral"] / normValue

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df

def create_dataframes_from_mresnova_json(data: dict)->dict:
    """
    Returns a dictionary of pandas dataframes for each technique in the data dictionary.
    """
    dataframes = {}
    for k, v in data.items():
        if k in ["H1_1D", "C13_1D", "HSQC", "HMBC", "COSY", "NOESY", "H1_pureshift"]:
            if v["type"].lower() == "2d":
                df = get_2D_dataframe_from_json(data, k)
                dataframes[k] = df
            elif v["type"].lower() == "1d":
                df = get_1d_dataframe_from_json(data, k)
                dataframes[k] = df
        elif k in ["smiles"]:
            dataframes["molecule"]  = pd.DataFrame([data["smiles"]], columns=["smiles"])

    return dataframes

def write_excel_file_from_mresnova_df_dict(df_frames:dict, excel_path: pathlib.Path)->bool:

    # check if path is valid
    if not excel_path.parent.exists():
        return False
    else:
        with pd.ExcelWriter(excel_path) as writer:
            # check if path is valid
            for k, df in df_frames.items():
                df.to_excel(writer, sheet_name=k)
        return True
