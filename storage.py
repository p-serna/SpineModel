# This is a class to store parameters of the experiments

class dataStorage(object): 
    def __init__(self,name='Data'):
        self.name = name
        
    def print(self):
        for k in self.__dict__.keys():
            print(k+" :",self.__dict__[k])
