import pandas as pd
import os
folder="FiguresData/"; 
afiles=os.listdir(folder)
files = []
for file in afiles:
    if file.split(".")[-1]=="csv":
        files.append(folder+file)

files.sort()
for file in files:
    temp = pd.read_csv(file)
    temp.to_excel(file.split(".")[0]+".xls")
