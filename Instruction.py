class Instruction:
    def __init__(self,Path):
        File=open(Path,'r')
        lines=File.read().split("\n")
        Data=[i.split(" ") for i in lines]
        for label in Data:
            if label[0]=="$":
                if len(label)>1:
                    if label[1]=='ENDFILE':
                        break
                continue
            elif label[0]=="TimeTot":
                self.TimeTot=int(label[1])
            elif label[0]=="seed":
                self.seed=int(label[1])
            elif label[0]=="B0":
                self.B0=float(label[1])
            elif label[0]=="Bf":
                self.Bf=float(label[1])
            elif label[0]=="SimulationType":
                self.SimulationType=label[1]
            elif label[0]=="Lx":
                self.Lx=int(label[1])
            elif label[0]=="Ly":
                self.Ly=int(label[1])
            elif label[0]=="mu":
                self.mu=float(label[1])
            elif label[0]=="MovingRate":
                self.MovingRate=float(label[1])
            elif label[0]=="TimeOut":
                self.TimeOut=int(label[1])
            elif label[0]=="J":
                self.J=float(label[1])
            elif label[0]=="TimeStat":
                self.TimeStat=int(label[1])
            elif label[0]=="Initial_System":
                self.Initial_System=label[1]
            elif label[0]=="NSimul":
                self.NSimul=int(label[1])
            elif label[0]=="k1":
                self.k1=float(label[1])
            elif label[0]=="k3":
                self.k3=float(label[1])
            elif label[0]=="k2":
                self.k2=float(label[1])
            elif label[0]=="Poisson":
                self.Poisson=float(label[1])
            elif label[0]=="eps":
                self.eps=float(label[1])
            elif label[0]=="Outputing":
                self.Outputing=label[1]
            elif label[0]=="kappa":
                self.kappa=label[1]
            elif label[0]=="Path":
                self.Path=label[1]
            elif label[0]=="Npart":
                self.Npart=label[1]
            else:
                print("data not recognize :")
                print(label)

