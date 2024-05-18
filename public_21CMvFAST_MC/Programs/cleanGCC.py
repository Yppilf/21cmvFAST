import os
import re

### This program is used to clean the output of the gcc compiler when using debug flags

files = [
    "ComputingTau_e.c",
    "create_dens_boxes_for_LC.c",
    "Createfcoll_ionisation_LC.c",
    "CreateFcollTable.c",
    "CreateSmoothedDensityBoxes.c",
    "drive_21cmMC_streamlined.c",
    "elec_interp.c",
    "filter.c",
    "folder_smooth.c",
    "heating_helper_progs.c",
    "init.c",
    "perturb_field.c",
    "smooth_output_box.c",
    "SplitMockObservation.c"
]

for file in files:
    compiledFiles = [f for f in os.listdir('.') if len(re.findall(f"{file}(\S+)", f)) != 0]
    for compiledFile in compiledFiles:
        os.remove(compiledFile)