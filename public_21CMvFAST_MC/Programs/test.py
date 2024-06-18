# import re

# class Test:
#     def __init__(self, walker_folder_location, IncludeLightCone, param_legend):
#         self.walker_folder_location = walker_folder_location
#         self.IncludeLightCone = IncludeLightCone
#         self.param_legend = param_legend
#         self.Fiducial_Params = {}
    
#     def read_walker(self):
#         # This is for dummy testing
#         params = [0]*100

#         walker_file = open(f"{self.walker_folder_location}Walker_1.000000_1.000000.txt","r")
#         Lines = walker_file.readlines()
#         self.parameter_number = 0

#         # FLAGS
#         match = re.match("^FLAGS\s+(\d)\s+(\d)\s+(\d)\s+(\d)\s+(\d)\s+(\d)\s+(\d)$", Lines[0].strip())
#         if match:
#             groups = match.groups()
#             GenerateNewICs = groups[0]
#             Subcell_RSDs = groups[1]
#             IONISATION_FCOLL_TABLE = groups[2]
#             UseFcollTable = groups[3]
#             PerformTsCalc = groups[4]
#             INHOMO_RECO = groups[5]
#             OutputGlobalAve = groups[6]
#         else:
#             raise Exception(f"Walker file: {walker_file} does not match expected pattern for FLAGS")

#         def read_float(string, saveKey, alwaysUseFiducial=False, useOtherParamKey=False, otherParamKey=None):
#             match = re.match(f"^{saveKey}\s+([0-9.]+)$", string.strip())
#             if match:
#                 groups = match.groups()
#                 if not alwaysUseFiducial:      
#                     if useOtherParamKey:
#                         paramKey = otherParamKey
#                     else:
#                         paramKey = saveKey     
#                     if self.param_legend[paramKey]:            
#                         params[self.parameter_number] = groups[0]
#                         self.parameter_number += 1
#                     else:
#                         self.Fiducial_Params[saveKey] = groups[0]   # Double because nested if cant be evaluated for some settings
#                 else:
#                     self.Fiducial_Params[saveKey] = groups[0]
#             else:
#                 raise Exception(f"Walker file: {walker_file} does not match expected pattern for {saveKey}")   

#         read_float(Lines[1], "ALPHA")
#         read_float(Lines[2], "ZETA")
#         read_float(Lines[3], "MFP")
#         read_float(Lines[4], "TVIR_MIN")
#         X_RAY_TVIR_MIN = params[self.parameter_number]
#         read_float(Lines[5], "L_X")
#         read_float(Lines[6], "NU_X_THRESH")
#         read_float(Lines[7], "NU_X_BAND_MAX", alwaysUseFiducial=True)
#         read_float(Lines[8], "NU_X_MAX", alwaysUseFiducial=True)
#         read_float(Lines[9], "X_RAY_SPEC_INDEX")
#         read_float(Lines[10], "X_RAY_TVIR_MIN", useOtherParamKey=True, otherParamKey="TVIR_MIN")
#         read_float(Lines[11], "X_RAY_TVIR_LB", alwaysUseFiducial=True)

#         read_float(Lines[12], "X_RAY_TVIR_UB", alwaysUseFiducial=True)
#         read_float(Lines[13], "F_STAR", alwaysUseFiducial=True)
#         read_float(Lines[14], "t_STAR", alwaysUseFiducial=True)
#         read_float(Lines[15], "N_RSD_STEPS", alwaysUseFiducial=True)
#         read_float(Lines[16], "LOS_direction", alwaysUseFiducial=True)

#         # Ts: I haven't added the CO_EVAL_Z here as I don't need it at the moment

#         parameter_number = self.parameter_number

# # walker_folder_location = "/scratch/s4950836/Walker/"
# walker_folder_location = "./"
# param_legend = {    # It doesnt matter what the values are here, just that they all exist
#     "ALPHA": True,
#     "ZETA": True,
#     "MFP": False,
#     "TVIR_MIN": False,
#     "L_X": True,
#     "NU_X_THRESH": False,
#     "X_RAY_SPEC_INDEX": True,
#     "TVIR_MIN": True,
# }
# test = Test(walker_folder_location, True, param_legend)
# test.read_walker()

import os

os.system("echo hi")