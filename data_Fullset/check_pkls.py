import pickle
import numpy as np

with open("Fullset.pkl", "rb") as f:
    fset = pickle.load(f)
with open("FullsetSR+Rcorr.pkl", "rb") as f:
    newfset = pickle.load(f)
with open("Fullset_shrnk_corrected.pkl", "rb") as f:
    fset_cr = pickle.load(f)

newfsetks = newfset.keys()
fsetks = fset.keys()
fset_crks = fset_cr.keys()

#print(fset_crks == fsetks)

print(len(list(fsetks)),len(list(newfsetks)))

for key in newfsetks:
    print(key,":", type(newfset.get(key)),type(fset.get(key)))

for key in newfsetks:
    if fset.get(key) is not None:
        print(key)
        print(((type(newfset.get(key).iloc[0]), type(fset_cr.get(key).iloc[0]))))


for key in newfsetks:
    if fset.get(key) is not None:
        print(key)
        print(np.column_stack((newfset.get(key).iloc[:4].values, fset_cr.get(key).iloc[:4].values)))

