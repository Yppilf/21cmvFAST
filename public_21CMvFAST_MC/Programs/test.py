import subprocess

command = "ls ./CosmoHammer_21CMMC/help"
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
process.wait()
print(process.returncode)