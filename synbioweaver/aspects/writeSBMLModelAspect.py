from synbioweaver.core import *
from synbioweaver.aspects.reactionDefinitions import *
import numpy, os, copy
from libsbml import *

# most of this code is hacked from teh developers tutorial
# http://sbml.org/Software/libSBML/docs/python-api/libsbml-python-creating-model.html

def check(value, message):
   """If 'value' is None, prints an error message constructed using
   'message' and then exits with status code 1.  If 'value' is an integer,
   it assumes it is a libSBML return status code.  If the code value is
   LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
   prints an error message constructed using 'message' along with text from
   libSBML explaining the meaning of the code, and exits with status code 1.
   """
   if value == None:
       raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
   elif type(value) is int:
       if value == LIBSBML_OPERATION_SUCCESS:
           return
       else:
           err_msg = 'Error encountered trying to ' + message + '.' \
                     + 'LibSBML returned error code ' + str(value) + ': "' \
                     + OperationReturnValue_toString(value).strip() + '"'
           raise SystemExit(err_msg)
   else:
       return

def addSpecies(model,sid):
    s = model.createSpecies()
    check(s,                                 'create species s1')
    check(s.setId(sid),                      'set species s1 id')
    check(s.setCompartment('c1'),            'set species s1 compartment')
    check(s.setConstant(False),              'set "constant" attribute on s1')
    check(s.setInitialAmount(0),             'set initial amount for s1')
    check(s.setSubstanceUnits('mole'),       'set substance units for s1')
    check(s.setBoundaryCondition(False),     'set "boundaryCondition" on s1')
    check(s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on s1')

def addParameter(model,pid):
    p = model.createParameter()
    check(p,                                  'create parameter k')
    check(p.setId(pid),                       'set parameter k id')
    check(p.setConstant(True),                'set parameter k "constant"')
    check(p.setValue(1),                      'set parameter k value')
    check(p.setUnits('per_minute'),           'set parameter k units')
 
def addReaction(model,rid,reactants,products,rateLaw):
    reac = model.createReaction()
    check(reac,                                 'create reaction')
    check(reac.setId(rid),                      'set reaction id')
    check(reac.setReversible(False),            'set reaction reversibility flag')
    check(reac.setFast(False),                  'set reaction "fast" attribute')

    for rt in reactants:
        species_refr = reac.createReactant()
        check(species_refr,                        'create reactant')
        check(species_refr.setSpecies(rt),         'assign reactant species')
        check(species_refr.setConstant(False),     'set "constant" on species ref 1')
        check(species_refr.setStoichiometry(1.0),  'set stochiometry')

    for pd in products:
        species_refp = reac.createProduct()
        check(species_refp,                        'create reactant')
        check(species_refp.setSpecies(pd),         'assign reactant species')
        check(species_refp.setConstant(False),     'set "constant" on species ref 1')
        check(species_refp.setStoichiometry(1.0),  'set stochiometry')

    math_ast = parseL3Formula(rateLaw)
    check(math_ast,                           'create AST for rate expression')
    kinetic_law = reac.createKineticLaw()
    check(kinetic_law,                        'create kinetic law')
    check(kinetic_law.setMath(math_ast),      'set math on kinetic law')


class WriteSBMLModel(Aspect):
    
   def mainAspect(self):
      self.addWeaverOutput(self.writeSBMLModel)

   def writeSBMLModel(self,weaverOutput):
      # We are expecting either a set of reactions or a context

      if getattr(weaverOutput, "getContext", None) != None:
         self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getContext()
      else:
         if getattr(weaverOutput, "getReactions", None) != None:
            self.nspecies, self.nreactions, self.species, self.reactions, self.stoichiometry_matrix, self.parameters = weaverOutput.getReactions()
         else:
            print "printReactionNetwork : Neither getContext() or getReactions() is available. Quitting"
            exit()

      return self.printSBMLModel()
        

   def printSBMLModel(self):
        
      try:
         document = SBMLDocument(3, 1)
      except ValueError:
         raise SystemExit('Could not create SBMLDocumention object')

      model = document.createModel()
      check(model,                              'create model')
      check(model.setTimeUnits("minute"),       'set model-wide time units')
      check(model.setExtentUnits("mole"),       'set model units of extent')
      check(model.setSubstanceUnits('mole'),    'set model substance units')
 
      per_minute = model.createUnitDefinition()
      check(per_minute,                         'create unit definition')
      check(per_minute.setId('per_minute'),     'set unit definition id')
      unit = per_minute.createUnit()
      check(unit,                               'create unit on per_minute')
      check(unit.setKind(UNIT_KIND_SECOND),     'set unit kind')
      check(unit.setExponent(-1),               'set unit exponent')
      check(unit.setScale(0),                   'set unit scale')
      check(unit.setMultiplier(60),             'set unit multiplier')
 
      # Create a compartment inside this model, and set the required
      # attributes for an SBML compartment in SBML Level 3.
 
      c1 = model.createCompartment()
      check(c1,                                 'create compartment')
      check(c1.setId('c1'),                     'set compartment id')
      check(c1.setConstant(True),               'set compartment "constant"')
      check(c1.setSize(1),                      'set compartment "size"')
      check(c1.setSpatialDimensions(3),         'set compartment dimensions')
      check(c1.setUnits('litre'),               'set compartment size units')

      # original test model
      #addSpecies(model, "s1")
      #addSpecies(model, "s2")
      #addParameter(model, "k")
      #addReaction(model, "r1", ["s1"], ["s2"], "k * s1 * c1")

      for sp in self.species:
         addSpecies(model, str(sp) )
         
      
      for rn in self.reactions:
         addReaction(model, "r", rn.reactants, rn.products, rn.rate)
         print rn.rate

      print self.parameters 
      for pt in self.parameters:
         addParameter(model, pt)
         
      return writeSBMLToString(document)
