import os
import nmrProblem

class simpleNMRresults:

    def __init__(self, nmrproblem):
        self.nmrproblem = nmrproblem

        self.problem_data = None
        self.xy3_data = None

        ok, msg = self.readNMRproblem()

        

    # read in nmrproblem and xy3 json files from the nmrproblem directory
    def readNMRproblem(self):
        
        fn_json = os.path.join(self.nmrproblem.problemDirectoryPath, f'{self.nmrproblem.problemDirectory}.json')
        fn_xy3_json = os.path.join(self.nmrproblem.problemDirectoryPath, 'xy3.json')

        if os.path.exists(fn_json):
            # read in the fn_json file
            with open(fn_json, 'r') as f:
                self.problem_data = json.load(f)
        else:
            print(f'Error: {fn_json} does not exist')
            return False, "Error: {fn_json} does not exist"

        if os.path.exists(fn_xy3_json):
            # read in the fn_xy3_json file
            with open(fn_xy3_json, 'r') as f:
                self.xy3_data = json.load(f)
        else:
            print(f'Error: {fn_xy3_json} does not exist')
            return False, "Error: {fn_xy3_json} does not exist"

        return True, "ok"



