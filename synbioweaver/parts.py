# list of standard parts

promoters = ["PTet"]

class dbPTet:
    def __init__(self):
        self.inputs = ['aTc']
        self

    def metaData(self):
        info ="""
        Transfer function obtained from Tamsir et. al. Nature (2010). 
        Assumes promoter has TetR bound and is then induced by aTc
        """

        return info

    def transferFunction(self):
        return "3000*350/(1 + 350 + 2*160*(1-aTc/(1.7 + aTc)) + 160*160*(1-aTc/(1.7 + aTc))*(1-aTc/(1.7 + aTc)) )"

    
