from synbioweaver.core import *
from synbioweaver.aspects.modelDefinitions import *
from synbioweaver.aspects.molecularReactions import *
import numpy, os
from copy import *

class MassActionKineticsRNA(Aspect, MolecularReactions):
    
    def mainAspect(self):
        MassActionKineticsRNA.builtReactions = False
        self.addWeaverOutput(self.getReactions)

        self.reactions = []
        self.species = []
        self.parameters = []

    def getReactions(self, weaverOutput):
        
        if getattr(weaverOutput, "buildPromoterMap", None) == None:
            sys.exit("MassActionKineticsRNA : buildPromoterMap() is unavailable. Quitting"  )

        if MassActionKineticsRNA.builtReactions == False:
            self.promoterMap = weaverOutput.buildPromoterMap()
          
            # get molecular reactions first
            mmol = MolecularReactions()
            mreactions, mspecies = mmol.getMolecules(weaverOutput)

             # get the unique namespaces and generate reactions
            nspaces = getNamespaces(self.promoterMap)
            for ns in nspaces:
                reactions = self.getReactionsMassActionRNA(ns)
                
                # assign mass action rates and parameters to all processes
                for r in reactions:
                    r.assignMassAction()
                    for k in r.param:
                        self.parameters.append( k )
                
                self.reactions = self.reactions + reactions

                # assign mass action rates to molecular reactions
                for r in mreactions:
                    # get either reactant or product
                    if len(r.reactants) != 0:
                        scope = r.reactants[0].scope
                    else:
                        scope = r.products[0].scope
        
                    if scope == ns:
                        r.assignMassAction()
                        for k in r.param:
                            self.parameters.append( k )
                        
            
            self.reactions = self.reactions + mreactions 
            self.species = self.species + mspecies

            self.newModel = Model( self.species, self.reactions, self.parameters )
           
            MassActionKineticsRNA.builtReactions = True

        return deepcopy( self.newModel )

    
    def getReactionsMassActionRNA(self, nspace):
        reactions = []

        for key in range( len(self.promoterMap) ):
            mapping = self.promoterMap[key]

            if mapping.getScope() != nspace:
                continue

            partname = mapping.getId()
            regulators = mapping.getRegulators()
            codings = mapping.getCodings()

            if len(regulators) == 0:
                # This is a const promoter
                # need to add:
                # pr -> mX + pr, mX -> X, mX -> 0, X -> 0 
                
                for p in codings:
                    self.species.append( Species(mapping.getScope(),str(p), "protein") )
                    self.species.append( Species(mapping.getScope(),"m"+str(p), "mRNA") )

                # production of mRNAs
                prods1 = [Species(mapping.getScope(),"m"+str(x)) for x in codings]
                reactions.append( Reaction([], prods1, "rnaTransc" ) )

                # translation
                for p in codings:
                    # translation
                    reactions.append( Reaction( [Species(mapping.getScope(),"m"+str(p))], [Species(mapping.getScope(),str(p))], "proteinTransl" ) )
                    # decay of mRNAs
                    reactions.append( Reaction( [Species(mapping.getScope(),"m"+str(p))], [], "rnaDeg") )
                    # decay of proteins
                    reactions.append( Reaction([Species(mapping.getScope(),str(p))], [], "proteinDeg") )

            else:
                for i in range(len(regulators)):
                    regulator = regulators[i]

                    complx = partname + '_' + str(regulator)

                    Spartname = Species( mapping.getScope(),partname, "promoter" )
                    Sregulator = Species( mapping.getScope(), str(regulator), "protein" )
                    Scomplx = Species( mapping.getScope(), complx, "bindingDNA" )
                
                    # the binding/unbinding reactions Pr + R <-> Pr_R
                    reactions.append( Reaction([Spartname, Sregulator], [Scomplx], "dnaBind" ) )
                    reactions.append( Reaction([Scomplx], [Spartname, Sregulator], "dnaUnbind") )
                    self.species.append( Spartname )
                    self.species.append( Sregulator )
                    self.species.append( Scomplx )

                    # if positive promoter then the bound complex expresses
                    if mapping.getPolarities()[i] == 1:
                        mprods = [Species(mapping.getScope(), "m"+str(x)) for x in codings]
                        mprods.append(Scomplx)

                        reactions.append( Reaction([Scomplx], mprods, "rnaTransc") )
                        for p in codings:
                            reactions.append( Reaction([Species(mapping.getScope(),p)], [], "proteinDeg") )
                            reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [], "rnaDeg") )
                            reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [Species(mapping.getScope(),p)], "proteinTransl") )
                            self.species.append( Species(mapping.getScope(),p, "protein") )
                            self.species.append( Species(mapping.getScope(),"m"+str(p), "mRNA") )

                    # if negative promoter then just the promoter expresses
                    else:
                        mprods = [Species(mapping.getScope(),"m"+str(x)) for x in codings]
                        mprods.append(Spartname)
                        reactions.append( Reaction([Spartname], mprods,"rnaTransc") )

                        for p in codings:
                            reactions.append( Reaction([Species(mapping.getScope(),p)], [], "proteinDeg") )
                            reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [], "rnaDeg") )
                            reactions.append( Reaction([Species(mapping.getScope(),"m"+str(p))], [Species(mapping.getScope(),p)], "proteinTransl") )
                            self.species.append( Species(mapping.getScope(),p, "protein") )
                            self.species.append( Species(mapping.getScope(),"m"+str(p), "mRNA") )

        return reactions
