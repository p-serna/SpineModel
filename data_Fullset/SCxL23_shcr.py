import json

with open("data_Fullset/SCxL23_shcr.json","r") as f:
    pars = json.load(f)

locals().update(pars)

#gtrG = 0.0058596073805616536 # 1 nS in average per synapse
#gtrN =  0.0034375 
#gtrA =  0.00315
#GABAtaus = (0.5,15.0)
#NMDAtaus = (0.5,17.0)
#dendDiam = 0.792
#rm = 12e3
#ra = 250
#factor2nddend = 75
#len2nddend = 70
#inhOutside = True
#factorspinesdend = 3.3
#dendsizeL0 = 10
#denddiam0 = 1.24

