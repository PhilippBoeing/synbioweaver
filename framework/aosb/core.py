from abc import ABCMeta, abstractmethod
import inspect
import warnings
import re
import types


class ExecutionNode(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        self.before = []
        self.after = []
        self.scope = None
        self.additionStack = None;

    def getBeforeNodes(self, filterType=object):
        result = []
        for executionNode in self.before:
            if isinstance(executionNode, filterType):
                result.append(executionNode)
        return result

    def getAfterNodes(self, filterType):
        result = []
        for executionNode in self.after:
            if isinstance(executionNode, filterType):
                result.append(executionNode)
        return result


# Parts are always used as objects
class Part(ExecutionNode):
    """Abstract Superclass for all Parts
    
    Parts are the atomistic instructions in the genetic circuit execution flow
    
    | *Attributes:*
    |     before : Part
    |         The part that immediately precedes this one in the execution flow
    |     after : Part
    |         The part that immediately follows this one in the execution flow
    |     moleculeConnection : class(Molecule)
    |         Each part may have a logical link to a molecule 
    |     namespace : Circuit or Aspect
    |         Each part belongs to a namespace, i.e. the circuit or aspect that has
    |         added it to the execution flow
    """

    __metaclass__ = ABCMeta
 #   namespace = None

    def __init__(self):
        super(Part,self).__init__()
        self.precompileMoleculesBefore = []
        self.precompileMoleculesAfter = []

    def setBeforePart(self, part):
        oldPart = self.getBeforePart()
        if oldPart != None:
            self.before.remove(oldPart)
        self.before.append(part)

    def setAfterPart(self, part):
        oldPart = self.getAfterPart()
        if oldPart != None:
            self.after.remove(oldPart)
        self.after.append(part)

    def getBeforePart(self):
        for executionNode in self.before:
            if isinstance(executionNode, Part):
                return executionNode
        return None

    def getAfterPart(self):
        for executionNode in self.after:
            if isinstance(executionNode, Part):
                return executionNode
        return None

    def weave(self, weaver):
        for beforeMol in self.precompileMoleculesBefore:
            mol = weaver.getMoleculeObject(self.scope, beforeMol)
            self.before.append(mol)
            mol.after.append(self)

        for afterMol in self.precompileMoleculesAfter:
            mol = weaver.getMoleculeObject(self.scope, afterMol,True)
            self.after.append(mol)
            mol.before.append(self)

    def __str__(self):
        return self.__class__.__name__

class MoleculeList(object):
    def __init__(self,first,second):
        moleculeList = [first,second]

    def __add__(self,other):
        moleculeList.append(other)

    def __rshift__(self, other):
        return MoleculeReaction(self,other)

class MoleculeReaction(object):
    def __init__(self,fromMol,toMol):
        self.fromMol = fromMol
        self.toMol = toMol

class Molecule(ExecutionNode):
    """Superclass for all molecules. 
    
    All molecules are used as classes, not objects. Any agents such as
    Proteins, Transcription Factors, etc. need to be subtypes of Molecules"""
   # namespace = None

    # before = []
    #after = [] 
    # for now these were removed, instead they are created when they are first added to the
    # moleculeList. This is, so that this static property belongs to the newly created class,
    # e.g. MoleculeA, GFP, etc... instead of the overall Molecule type
    # eventually, Molecules should be treated as objects, similarly to Parts. 
    # Although I don't yet like the naming convention for parts....

    def __str__(self):
        return self.__class__.__name__

    def __add__(self, other):
        return MoleculeList(self,other)

    def __rshift__(self, other):
        return MoleculeReaction(self,other)


class Protein(Molecule):
    """Class for Molecules that are Proteins"""
    pass


class Promoter(Part):
    """Class for Parts that are Promoters"""
    pass


class ConstitutivePromoter(Promoter):
    """Class for Promoters that are Constitutive Promoters"""
    pass


def checkIfTypeReturnInstance(possibleType):
    """if the parameter is a type, try to return an instance
    
    | *Args:*
    |     possibleType - a parameter which may be an instance or a type
    
    | *Raises:*
    |     PartInitializationError - if the parameter is a type which can not be constructed
    
    | *Returns:*
    |     An instance of the type of possibleType
    
    """
    if isinstance(possibleType, type):
        try:
            return possibleType()
        except:
            raise PartInitializationError(str(possibleType) + " can not be initialized without Parameters.")

    return possibleType


def checkAndSetMolecule(molecule):
    """checks if parameter is a class of type Molecule and returns it
    
    | *Args:*
    |     molecule - A potential Molecule
    
    | *Raises:*
    |     MoleculeValue - If molecule is not a class of type Molecule
    
    | *Returns:*
    |     A class of type Molecule
    """

    if inspect.isclass(molecule) and issubclass(molecule, Molecule):
        return molecule
    else:
        raise MoleculeValueError("molecule must be a class of (sub)type Molecule")



class NegativePromoter(Promoter):
    """Class for Promoters that are Negative Promoters"""

    def __init__(self, regulatedBy):
        """NegativePromoter Constructor
        
        | *Args:*
        |     regulatedBy: sets moleculeConnection
        """

        super(NegativePromoter, self).__init__()
        self.precompileMoleculesBefore.append( checkAndSetMolecule(regulatedBy) )

    def getRegulatedBy(self):
        result = self.getBeforeNodes(Molecule)
        if len(result) > 0:
            return result

        #else: not combiled yet
        return self.precompileMoleculesBefore



    def __str__(self):
        tmpstr =  super(NegativePromoter, self).__str__() + '(regulatedBy = '
        tmpstr +=  'm'.join(str(self.getBeforeNodes(Molecule)[i]) for i in range(len(self.getBeforeNodes(Molecule))))
        tmpstr += ')'
        return tmpstr


class PositivePromoter(Promoter):
    """Class for Promoters that are Positive Promoters"""

    def __init__(self, regulatedBy):
        """PositivePromoter Constructor
        
        | *Args:*
        |     regulatedBy: sets moleculeConnection
        """

        super(PositivePromoter, self).__init__()
        self.precompileMoleculesBefore.append( checkAndSetMolecule(regulatedBy) )

# todo : is this necessary?
#    def __del__(self):
 #       self.regulatedBy.after.remove(self)

    def getRegulatedBy(self):
        result = self.getBeforeNodes(Molecule)
        if len(result) > 0:
            return result

        #else: not combiled yet
        return self.precompileMoleculesBefore

    def __str__(self):
        tmpstr =  super(PositivePromoter, self).__str__() + '(regulatedBy = '
        tmpstr +=  ','.join(str(self.getBeforeNodes(Molecule)[i]) for i in range(len(self.getBeforeNodes(Molecule))))
        tmpstr += ')'
        return tmpstr





class CodingRegion(Part):
    """Class for parts that are Coding Regions"""

    def __init__(self, codesFor):
        """Constructor for a Coding Region Part
        
        | *Args:* 
        |     codesFor: sets moleculeConnection
        """
        super(CodingRegion, self).__init__()
        self.precompileMoleculesAfter.append( checkAndSetMolecule(codesFor) )



 # todo : is this necessary?
#    def __del__(self):
 #       self.codesFor.before.remove(self)


    def getCodingFor(self):
        result = self.getAfterNodes(Molecule)
        if len(result) > 0:
            return result

        #else: not combiled yet
        return self.precompileMoleculesAfter

    def __str__(self):
        tmpstr =  super(CodingRegion, self).__str__() + '(codesFor = '
        tmpstr +=  'm'.join(str(self.getAfterNodes(Molecule)[i]) for i in range(len(self.getAfterNodes(Molecule))))
        tmpstr += ')'
        return tmpstr


class RBS(Part):
    """Class for parts that are Ribosome BindingSites"""
    pass


class Terminator(Part):
    """Class for parts that are Terminators"""
    pass


class Circuit(Part):
    """Abstract class for a genetic part circuit.
    
    | Will generally be used to represent a design's core concerns
    | Additionally, a circuit can be used to represent composite parts
    
    | *Attributes:*
    |     weaver: Weaver
    |         The weaver which compiles the circuit
    """

    __metaclass__ = ABCMeta
    weaver = None

    @abstractmethod
    def mainCircuit(self):
        """Entry point for a circuit, analogous to "main" in a program
       
        | Needs to be implemented by any sub class.
        | mainCircuit will be called by the AOSB Weaver."""
        pass

    def importMolecule(self, molecule):
        self.weaver.importMolecule(self,molecule)

    def exportMolecule(self, molecule):
        self.weaver.exportMolecule(self,molecule)

    def createMolecule(self, molecule):
        self.weaver.createMolecule(self,molecule)

    def addCircuit(self, circuit):
        self.weaver.addCircuit(self,circuit)

    def addPart(self, part):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver
        
        | *Args:*
        |     part: The part to be added to the circuit
        
        """

        self.weaver.addPart(self, part)

    def reactionFrom(self, *molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """

        return self.weaver.reactionFrom(self, molecules)

    def reactionTo(self, *molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """

        return self.weaver.reactionTo(self, molecules)

    def setWeaver(self, weaver):
        """Internal - Should not be used outside of the framework.
        
        Sets the Weaver Object of this Circuit
        
        | *Args:*
        |     weaver: A weaver object that will be used
            """

        self.weaver = weaver


def declareNewPart(classname, parent=Part, moleculesBefore=[], moleculesAfter = []):
    ''' Returns a new Part type and exports it to the caller's namespace
    
    | *Args*
    |     classname : string
    |         The name for the new type
    |     parent : Part
    |         super class for the new type
    |     moleculeConnection : Molecule
    |         optional, if the new Part type should be connected to a Molecule type
    
    | *Returns*
    |     The new Part type
    
    | *Raises*
    |     PartValueError -If the parent is not a Part
    |     InvalidSymbolNameError - If classname is not a valid name for a symbol
    
    | *Warnings*
    |     SymbolExistsWarning - If classname already exists in the namespace'''


    if not issubclass(parent,Part):
        raise PartValueError("parent must be of type Part")   
    
    validnameregex = re.compile('[a-zA-Z_][a-zA-Z0-9_]*')
        
    if not isinstance(classname,str) or not validnameregex.match(classname) or not validnameregex.match(classname).span() == (0,len(classname)):
        raise InvalidSymbolNameError('name is not a valid symbold name')

    if not isinstance(moleculesBefore,list):
        raise ValueError("moleculesBefore must be of type list")

    if not isinstance(moleculesAfter,list):
        raise ValueError("moleculesAfter must be of type list")

    for mol in moleculesBefore:
        checkAndSetMolecule(mol)

    for mol in moleculesAfter:
        checkAndSetMolecule(mol)

    basestuple = parent,
    result = None
    currentframe = inspect.currentframe()
    # Create warning if name already exists 
    if currentframe.f_back.f_globals.has_key(classname):
        line = currentframe.f_back.f_lineno
        warnings.warn("Line "+str(line)+": Part "+classname+" already defined. Existing definition will be used",SymbolExistsWarning,2) # 2 = one stack level above this
        result = currentframe.f_back.f_globals.get(classname)
    else:
        result = type(classname,basestuple,{})
        def newTypeInit(self,recursion = 0):
            if(recursion == 0):
                Part.__init__(self)
            try:
                realself = self.__class__
                for i in range(0,recursion):
                    realself = realself.__bases__[0]

                super(realself,self).__init__(recursion+1)
            except:
                pass

            self.precompileMoleculesBefore += moleculesBefore
            self.precompileMoleculesAfter += moleculesAfter

        #this means that parent does not take a moleculeConnection parameter,
        #but the new subtype should have one.
        result = type(classname,basestuple,{'__init__':newTypeInit})
        currentframe.f_back.f_globals[classname] = result
    return result



def declareNewMolecule(classname, *parents):
    """Returns a new Molecule type and exports it to the caller's namespace
    
    | *Args:*
    |     classname: The name for the new type
    |     \*parents: 0 or more Molecule super classes
    
    | *Returns:*
    |     The new Part type
    
    | *Raises:*
    |     MoleculeValueError: If any parent is not a Molecule
    |     InvalidSymbolNameError: If classname is not a valid name for a symbol
    
    | *Warnings:*
    |     SymbolExistsWarning: If classname already exists in the namespace
    """

    validnameregex = re.compile('[a-zA-Z_][a-zA-Z0-9_]*')

    if not isinstance(classname, str) or not validnameregex.match(classname) or not validnameregex.match(
            classname).span() == (0, len(classname)):
        raise InvalidSymbolNameError('name is not a valid symbold name')

    if len(parents) == 0:
        parents = Molecule,
    else:
        for parent in parents:
            if not issubclass(parent, Molecule):
                raise MoleculeValueError("All parents must be of type Molecule")
    currentframe = inspect.currentframe()

    result = None

    # Check if classname already exists 
    if currentframe.f_back.f_globals.has_key(classname):
        line = currentframe.f_back.f_lineno
        warnings.warn(
            "Line " + str(line) + ": Molecule " + classname + " already defined. Existing definition will be used",
            SymbolExistsWarning, 2)  # 2 = one stack level above this
        result = currentframe.f_back.f_globals.get(classname)
    else:
        result = type(classname, parents, {})
        currentframe.f_back.f_globals[classname] = result
    return result


class PointCutExpressionNode(object):
    """Abstract superclass for all Nodes in a Point Cut Expression Tree
    
    | A complex expression for a Point Cut, using operators such as & (and),
    | \|  (or) or % (concatenation) is represented as a tree of nodes

    """

    __metaclass__ = ABCMeta

    def __mod__(self, other):
        """% - concatenates two PointCutExpressionNodes
        
        | *Args:*
        |     other - the node on the right of self
        
        | *Returns:*
        |     a new PointCutExpressionConcatenate node
        """

        return PointCutExpressionConcatenate(self, other)

    def __and__(self, other):
        """& - boolean and evaluation of two PointCutExpressionNodes
        
        | *Args:*
        |     other - the node on the right of self
        
        | *Returns:*
        |     a new PointCutExpressionNodeAnd node
        """

        return PointCutExpressionAnd(self, other)

    def __or__(self, other):
        """|  - boolean or evaluation of two PointCutExpressionNodes
        
        | *Args:*
        |     other - the node on the right of self
        
        | *Returns:*
        |     a new PointCutExpressionOr node
        """

        return PointCutExpressionOr(self, other)

    def numberOfMatchingParts(self, part):
        """returns the number of parts the expression matches
        | If an expression uses concatenation, then it might match the current part and a number of preceding parts.
        
        | *Args:*
        |     part - The part at which matching starts
        
        | *Returns:*
        |     integer - how many parts were matched
        """

        if self.match(part):
            return 1
        else:
            return 0

    @abstractmethod
    def match(self, part):
        """Does a part match this (sub)-expression? An abstract method, must be implemented by each child
        
        | *Args:*
        |     part - The part to be matched
        
        | *Returns:*
        |     boolean - Whether or not the part was matched
        """

        pass


class PointCutExpressionOperator(PointCutExpressionNode):
    """Abstract superclass of an PointCutExpression Node which is an operator
    
    | *Attributes:*
    |     left - The first child of the operator
    |     right - The second child of the operator
    """

    __metaclass__ = ABCMeta

    left = None
    right = None

    def __init__(self, left, right):
        """Constructs a new PointCutExpresionOperator
        
        | *Args:*
        |     left - Sets the left child node
        |     right - Sets the right child node
        
        | *Raises:*
        |     InvalidPointCutExpressionError - 
        |         If either child is not an instance of PointCutExpressionNode
        """

        if isinstance(left, PointCutExpressionNode) and isinstance(right, PointCutExpressionNode):
            self.left = left
            self.right = right
        else:
            raise InvalidPointCutExpressionError("Invalid type used in a PointCut formula.")

    def expressionUses(self, nodeType):
        """Confirms if a certain type of node is used in the expression
        
        | *Args:*
        |     nodeType: The type of the node whose existence is to be confirmed
        
        | *Returns:*
        |     boolean - Whether or not the node exists in the formula
        """

        if isinstance(self, nodeType) or isinstance(self.left, nodeType) or isinstance(self.right, nodeType):
            return True
        else:
            leftresult = False
            rightresult = False
            try:
                leftresult = self.left.expressionUses(nodeType)
            except:
                pass
            try:
                rightresult = self.right.expressionUses(nodeType)
            except:
                pass
            return leftresult or rightresult


class PointCutExpressionNot(PointCutExpressionOperator):
    """A PointCutExpressionOperator which is a Not
    
    A special case, only uses one child, acts as the inverse operator
    """

    def __init__(self, pointcutexpression):
        """Constructs a new PointCutExprresionNot
        
        | *Args:*
        |     pointcutexpression - The expression to be negated
        """

        self.right = pointcutexpression

    def match(self, part, within = None):
        """see PointCutExpressionNode definition"""

        return not self.right.match(part,within)


class PointCutExpressionOr(PointCutExpressionOperator):
    """A PointCutExpressionOperator which is an Or"""

    def match(self, part, within = None):
        """see PointCutExpressionNode definition"""

        return self.left.match(part,within) or self.right.match(part,within)


class PointCutExpressionAnd(PointCutExpressionOperator):
    """A PointCutExpressionOperator which is an And"""

    def match(self, part):
        """see PointCutExpressionNode definition"""

        return self.left.match(part) and self.right.match(part)


class PointCutExpressionConcatenate(PointCutExpressionOperator):
    """A PointCutExpressionOperator which is an Concatenation"""

    def numberOfMatchingParts(self, part):
        """see PointCutExpressionNode definition
        
        This child overrides it, since a concatenation operator is
        a node in the expression at which more than one part can be
        matched. 
        """

        if self.match(part):
            return self.left.numberOfMatchingParts(part.getBeforePart()) + self.right.numberOfMatchingParts(part)
        else:
            return 0


    def match(self, part):
        """see PointCutExpressionNode definition"""

        try:
            return self.left.match(part.getBeforePart()) and self.right.match(part)
        except:
            return False


def checkBaseClassesMatch(bases, typename):
    """Recursively check if the name of the types in bases (or parents) are equal to typename
    
    | *Args:*
    |     bases - A tuple of types
    |     typename : str - A name of a type
    
    | *Returns:*
    |     True if any of the names of types in bases or any of their parent bases equals typename,
    |         False otherwise
    """

    for baseClass in bases:
        if (typename == baseClass.__name__):
            return True
        else:
            if checkBaseClassesMatch(baseClass.__bases__, typename):
                return True


class PartSignatureElement(object):
    """A building block of a Part Signature
    
    | *Attributes:*
    |     qualifier : ANY / SUBCLASS / CLASSONLY
    |         Whether the element should match precisely, all subclasses or uses a wildcard
    |     element : str - The string of the element (without a qualifier)
    |     inverse : boolean - Whether the element has been negated
    """

    ANY = 1
    SUBCLASS = 2
    CLASSONLY = 3


    def __init__(self, signature):
        """Constructs a new PartSignatureElement 
        | Analyzes signature and sets internal attributes
        
        | *Args:*
        |     signature - the part of the Part signature for this element
        """

        self.inverse = False
        self.qualifier = None
        self.element = ''

        if signature[0] == '!':
            self.inverse = True
            signature = signature[1:]

        if signature[-1] == '*':
            self.qualifier = PartSignatureElement.ANY
            self.element = signature[:-1]
        elif signature[-1] == '+':
            self.qualifier = PartSignatureElement.SUBCLASS
            self.element = signature[:-1]
        else:
            self.qualifier = PartSignatureElement.CLASSONLY
            self.element = signature

    def __str__(self):
        qualifier = ''
        if self.qualifier == self.ANY:
            qualifier = '*'
        elif self.qualifier == self.SUBCLASS:
            qualifier = '+'
        exclamation = ''
        if self.inverse:
            exclamation = '!'
        return exclamation + self.element + qualifier

    def match(self, obj):
        """see PointCutExpressionNode definition"""

        if self.inverse:
            return not self.__match(obj)
        else:
            return self.__match(obj)

    def __match(self, obj):
        """internal match method, matching without inverse"""


        objectName = obj.__class__.__name__
        if(isinstance(obj,str)):
            objectName = obj
        if(inspect.isclass(obj)):
            objectName = obj.__name__

        if ((self.qualifier == self.CLASSONLY) or (self.qualifier == self.SUBCLASS)) and (objectName == self.element):
            return True

        if (self.qualifier == self.ANY) and (objectName.startswith(self.element)):
            return True

        if (self.qualifier == self.SUBCLASS):
            if(inspect.isclass(obj)):
                return checkBaseClassesMatch(obj.__bases__, self.element)
            else:
                return checkBaseClassesMatch(obj.__class__.__bases__, self.element)

        return False


class MoleculeSignatureElement(PartSignatureElement):
    """A PartSignatureElement that is a Molecule Signature
    | Overloads some methods, since MoleculeSignature is used by the user for
    | Molecule Type Advice - unlike PartSignatureElement, which is internal
    """

    def __init__(self, signature):
        """Constructs a new MoleculeSignature
        
        | *Args:*
        |     signature : str - The string signature to be cast to a MoleculeSignature
        
        | *Raises:*
        |     InvalidSignatureError - If signature does not adhere to the format
        """

        signatureRE = re.compile('!?([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?|\*)')

        try:
            if not signatureRE.match(signature).span() == (0, len(signature)):
                raise InvalidSignatureError()
        except:
            raise InvalidSignatureError()

        return super(MoleculeSignatureElement, self).__init__(signature)

    def match(self, obj):
        """see PointCutExpressionNode definition"""

        if self.qualifier == self.ANY and obj == None:
            if self.inverse:
                return False
            return True
        return super(MoleculeSignatureElement, self).match(obj)

class MoleculeSignature():
    """A MoleculeSignature, used for Type Advice of Molecules

    | *Attributes:*
    |     namespace : [PartSignatureElement] - List of signature parts of the signature
    |     molecule : MoleculeSignatureElement - The "part" part of the signature
    """


    def __init__(self, signature):
        """Constructs a new PartSignature
        | Analyzes signature and sets internal attributes

        | *Args:*
        |     signature : str - The string signature to be cast to a PartSignature

        | *Raises:*
        |     InvalidSignatureError - If signature does not adhere to the format
        """

        self.namespace = []
        self.molecule = None

        # should be of format: Circuit.Part(Molecule)" (Molecule) is optional#
        typeErrorMessage = "signature must be of type String and adhere to PointCut / PartSignature Format"

        if not isinstance(signature, str):
            raise InvalidSignatureError(typeErrorMessage)


        # using a regular expression to make sure the signature is of a valid format
        signatureRE = re.compile(
            '(!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)\.)+!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)')


        # regular expression has to match entire length of string
        try:
            if not signatureRE.match(signature).span() == (0, len(signature)):
                raise InvalidSignatureError(typeErrorMessage)
        except:
            raise InvalidSignatureError(typeErrorMessage)

    #split signature into components
        signaturePartRE = re.compile('!?[a-zA-Z_][a-zA-Z0-9_]*[\+\*]?|\*')

        splitSignature = signaturePartRE.findall(signature)

        numOfCircuitElements = len(splitSignature)-1;


        self.molecule = MoleculeSignatureElement(splitSignature[numOfCircuitElements])
       # #Circuit Part
       # self.namespace = PartSignatureElement(splitSignature[0])

        #Part ... Part

        for i in range(0,numOfCircuitElements):
            self.namespace.append(PartSignatureElement(splitSignature[i]))

    def __str__(self):

        firsthalf =  '.'.join(str(self.namespace[i]) for i in range(len(self.namespace)))

        firsthalf = firsthalf + '.' + str(self.molecule)


    def matchNamespaces(self, namespace, scope):
        for name in reversed(namespace):
            if name.match(scope.circuitName):
                scope = scope.parentCircuit
            else:
                return False

        # todo edge case, i.e. if need to check * is last, or no match?
        return True

    def match(self, molecule):
        """see PointCutExpressionNode definition"""

        if not isinstance(molecule, Molecule):
            raise MoleculeValueError("Molecule to match must be instance of type Molecule")
        result = self.molecule.match(molecule) and self.matchNamespaces(self.namespace,molecule.scope)

        return result


class NotOnStack(PointCutExpressionNode):
    def __init__(self, signature):
        self.name = signature
        #todo all the error checking...

    def match(self, part):
        """see PointCutExpressionNode definition"""

        #todo matching here should also support + and *

        for withinObject in part.additionStack:
            if withinObject.__class__.__name__ == self.name:
                return False

        return True

class OnStack(PointCutExpressionNode):
    def __init__(self, signature):
        self.name = signature
        #todo all the error checking...

    def match(self, part):
        """see PointCutExpressionNode definition"""

        for withinObject in part.additionStack:
            if withinObject.__class__.__name__ == self.name:
                return True

        return False


class PartSignature(PointCutExpressionNode):
    """A PartSignature, used in PointCut expressions or directly for Type Advice
    
    | *Attributes:*
    |     namespace : [PartSignatureElement] - List of signature parts of the signature
    |     part : PartSignatureElement - The "part" part of the signature
    |     molcule : PartSignatureElement - The molecule part of the signature
    |     nomolecule : Boolean - 
    |         If the PartSignature explicitly should not match parts with molecules
    
    """


    def __init__(self, signature):
        """Constructs a new PartSignature
        | Analyzes signature and sets internal attributes
        
        | *Args:*
        |     signature : str - The string signature to be cast to a PartSignature
        
        | *Raises:*
        |     InvalidSignatureError - If signature does not adhere to the format
        """

        self.namespace = []
        self.part = None
        self.molecule = None
        self.nomolecule = False

        # should be of format: Circuit.Part(Molecule)" (Molecule) is optional#
        typeErrorMessage = "signature must be of type String and adhere to PointCut / PartSignature Format"

        if not isinstance(signature, str):
            raise InvalidSignatureError(typeErrorMessage)


        # using a regular expression to make sure the signature is of a valid format
        signatureRE = re.compile(
            '(!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)\.)+!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)(\(!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)?\))?')


        # regular expression has to match entire length of string
        try:
            if not signatureRE.match(signature).span() == (0, len(signature)):
                raise InvalidSignatureError(typeErrorMessage)
        except:
            raise InvalidSignatureError(typeErrorMessage)

    #split signature into components
        signaturePartRE = re.compile('!?[a-zA-Z_][a-zA-Z0-9_]*[\+\*]?|\*')

        splitSignature = signaturePartRE.findall(signature)

        numOfCircuitElements = len(splitSignature)-1;

        signatureMoleculePartRE = re.compile('\(!?(([a-zA-Z_][a-zA-Z0-9_]*[\+\*]?)|\*)?\)')



        moleculePartMatch = signatureMoleculePartRE.search(signature)
        if moleculePartMatch != None:
            moleculePart = moleculePartMatch.group(0)
            if len(moleculePart) == 2:
                self.nomolecule = True
            else:
                numOfCircuitElements -= 1
                self.molecule = MoleculeSignatureElement(moleculePart[1:-1])
                self.nomolecule = False

        else:
            self.nomolecule = False
            self.molecule = MoleculeSignatureElement('*')



        self.part = PartSignatureElement(splitSignature[numOfCircuitElements])
       # #Circuit Part
       # self.namespace = PartSignatureElement(splitSignature[0])

        #Part ... Part

        for i in range(0,numOfCircuitElements):
            self.namespace.append(PartSignatureElement(splitSignature[i]))





    def __str__(self):

        firsthalf =  '.'.join(str(self.namespace[i]) for i in range(len(self.namespace)))



        firsthalf = firsthalf + '.' + str(self.part)
        if self.nomolecule:
            return firsthalf + '()'
        else:
            return firsthalf + '(' + str(self.molecule) + ')'

    def matchNamespaces(self, namespace, scope):
        for name in reversed(namespace):
            if name.match(scope.circuitName):
                scope = scope.parentCircuit
            else:
                return False

        # todo edge case, i.e. if need to check * is last, or no match?
        return True

    def match(self, part):
        """see PointCutExpressionNode definition"""

        if not isinstance(part, Part):
            raise PartValueError("Part to match must be instance of type Part")
        result = self.part.match(part) and self.matchNamespaces(self.namespace,part.scope)

        if not self.nomolecule:
            # I.E. THERE IS A MOLECULE
            # this needs to be refined / TODO
            if len(part.precompileMoleculesAfter) == 0 and len(part.precompileMoleculesBefore) == 0:
                if self.molecule.qualifier == self.molecule.ANY and self.molecule.element == '':
                    result = result and True
                else:
                    result = result and False

            for mol in part.precompileMoleculesAfter:
                result = result and self.molecule.match(mol)

            for mol in part.precompileMoleculesBefore:
                result = result and self.molecule.match(mol)

        if self.nomolecule:
            if len(part.getBeforeNodes(Molecule)) != 0 or len(part.getAfterNodes(Molecule)) != 0 or\
                            len(part.precompileMoleculesAfter) != 0 or len(part.precompileMoleculesBefore) != 0:
                result = result and False

            #todo this means we could also do more than one molecule, e.g. Promoter(AB || A) or bool and...

        return result


class PointCut(object):
    """A PointCut to select Join Points in the genetic parts execution flow
    
    | *Attributes:*
    |     operator : BEFORE / AFTER / REPLACE 
    |         - The operator for this PointCut
    |     signature : PartSignature
    """
    BEFORE = 11
    AFTER = 22
    REPLACE = 33


    def __init__(self, signature, operator):
        """Construct a new PointCut
        
        | *Args:* 
        |     signature - A PointCutSignature or a string that will be cast to PointCutSignature
        |     operator - The operator
        """

        self.operator = None
        self.signature = None

        if isinstance(signature, PointCutExpressionNode):
            self.signature = signature
        else:
            self.signature = PartSignature(signature)
        self.checkAndSetOperator(operator)

    def __str__(self):
        operator = ''
        if self.operator == self.BEFORE:
            operator = 'BEFORE'
        elif self.operator == self.AFTER:
            operator = 'AFTER'
        else:
            operator = 'REPLACE'

        return 'PointCut(\'' + str(self.signature) + '\',' + operator + ')'

    def match(self, part):
        """see PointCutExpressionNode definition"""

        return self.signature.match(part)

    def checkAndSetOperator(self, operator):
        """Set the internal operator attribute, if the parameter is a valid operator
        
        | *Args:*
        |     operator - The operator to be confirmed
        
        | *Raises:*
        |     InvalidPointCutOperatorError - if the operator is invalid
        |         (can be dependent on the signature)
        """

        if (operator == self.BEFORE) or (operator == self.AFTER) or (operator == self.REPLACE):
            self.operator = operator
            if (operator == self.BEFORE and isinstance(self.signature,
                                                       PointCutExpressionOperator) and self.signature.expressionUses(
                    PointCutExpressionConcatenate)):
                raise InvalidPointCutOperatorError(
                    'PointCut Operator can not be BEFORE if PartSignature uses concatenation.')
        else:
            raise InvalidPointCutOperatorError("operator must be of type PointCut.BEFORE / AFTER / REPLACE")


class PointCutContext(object):
    """Container for context at a PointCut
    
    | *Attributes:*
    |     within - A stack of the circuits / aspects that within which the PointCut was matched
    |     part - The part matched by the PointCut"""

    def __init__(self, within, part):
        self.within = within
        self.part = part

    def isWithin(self, obj):
        """Checks if the parameter (a circuit or aspect) is on the within stack
        | *Args:*
        |     obj : Circuit / Aspect - The object that is to be found on the stack
        
        | *Returns:*
        |     Boolean - true if obj is on the within stack
        """
        return self.__isWithinRecursive(obj, len(self.within))

    def __isWithinRecursive(self, obj, index):
        if index > 0:
            return (obj == self.within[index - 1]) or self.__isWithinRecursive(obj, index - 1)
        return False


class Advice(object):
    """Container for Advice
    
    | *Attributes:*
    |     precedence : int - High precedence advice have execution priority over low precedence
    |         - MINPRECEDENCE <= precedence <= MAXPRECEDENCE is the valid range
    |     pointcut : PointCut
    |     adviceMethod : method - The method to be executed at the advice
    |         - must have 2 parameters: self and PointCutContext
    """

    MINPRECEDENCE = 0
    MAXPRECEDENCE = 100


    def __init__(self, pointcut, adviceMethod, precedence=MINPRECEDENCE):
        """Constructs a new Advice object, setting internal state
        
        | *Args:*
        |     pointcut - The PointCut to be set
        |     adviceMethod - The adviceMethod to be set
        |     precedence - The precedence to be set
        
        | *Raises:*
        |     InvalidPointCutError
        |     InvalidAdviceMethodError - If argument is not a method or takes wrong
        |         number of parameters (must take self and PointCutContext object)
        |     PrecedenceOutOfRangeError - If precedence > MAXPRECEDENCE or
        |                                     precedence < MINPRECEDENCE
        """

        self.pointcut = None
        self.adviceMethod = None
        self.precedence = self.MINPRECEDENCE

        if not isinstance(pointcut, PointCut):
            raise InvalidPointCutError("pointcut must be of type PointCut")
        self.pointcut = pointcut

        if not inspect.ismethod(adviceMethod) or len(inspect.getargspec(adviceMethod)[0]) != 2:
            raise InvalidAdviceMethodError("adviceMethod must be a method with 2 parameters (self, PointCutContext)")
        self.adviceMethod = adviceMethod

        if precedence < (self.MINPRECEDENCE) or (precedence > self.MAXPRECEDENCE):
            raise PrecedenceOutOfRangeError(
                "Advice Precedence must be between " + str(self.MINPRECEDENCE) + " and " + str(self.MAXPRECEDENCE))
        self.precedence = precedence


class TypeAdvice(object):
    """Container for TypeAdvice
    
    | *Attributes:*
    |     signature : PartSignature or MoleculeSignature
    |     typeaddition : method or attribute to be added to the type
    |     name : The name the new typeaddition should have in the new type
    |     aspect : Aspect which declares this TypeAdvice
    """


    def __init__(self, signature, typeaddition, name, aspect):
        """Constructs a new TypeAdvice
        
        | *Args:*
        |     signature - PartSignature / MoleculeSignature to be set
        |     typeaddition - attribute to be set
        |     name : str - attribute to be set
        |     aspect : Aspect - attribute to be set
        
        | *Raises:*
        |     InvalidSignatureError - 
        |         If signature not instance of PartSignature or MoleculeSignature
        |     InvalidSymbolNameError -
        |         If name is not a valid name for a symbol
        """

        self.signature = None
        self.typeaddition = None
        self.aspect = None
        self.name = ''

        if isinstance(signature, PartSignature) or isinstance(signature, MoleculeSignature):
            self.signature = signature;
        else:
            raise InvalidSignatureError(
                'Signature parameter for TypeAdvice must be of type PartSignature or MoleculeSignature')
        self.typeaddition = typeaddition;
        self.aspect = aspect

        validnameregex = re.compile('[a-zA-Z_][a-zA-Z0-9_]*')

        if not isinstance(name, str) or not validnameregex.match(name) or not validnameregex.match(name).span() == (
                0, len(name)):
            raise InvalidSymbolNameError('name is not a valid symbold name')
        self.name = name

    def isPartAdvice(self):
        """Returns True if TypeAdvice is for Part, False otherwise"""

        return isinstance(self.signature, PartSignature)

    def isMoleculeAdvice(self):
        """Returns True if TypeAdvice is for Molecule, False otherwise"""

        return isinstance(self.signature, MoleculeSignature)


class Aspect(object):
    """Abstract class for an Aspect
    | Will generally be used to represent a design's cross-cutting concerns
    
    | *Notes:*
    |     Any child must implement mainAspect method
    
    | *Attributes*
    |     weaver: Weaver
    |         The weaver which compiles the aspect
    |     adviceList - A list of all advice this aspect declares
    |     typeAdviceList - A list of all type advice this aspect declares
    |     weaverOutputList - A list of all type advice for the WeaverOutput
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        self.adviceList = []
        self.typeAdviceList = []
        self.weaverOutputList = []
        self.weaver = None

    def importMolecule(self, molecule):
        self.weaver.importMolecule(self,molecule)

    def exportMolecule(self, molecule):
        self.weaver.exportMolecule(self,molecule)

    def createMolecule(self, molecule):
        self.weaver.createMolecule(self,molecule)

    def addCircuit(self, circuit):
        self.weaver.addCircuit(self,circuit)

    def addPart(self, part):
        """Used to add parts to the design by passing them on to the AOSB Weaver
        
        | *Args:*
        |     part: The part to be added to the circuit
        
        """
        self.weaver.addPart(self, part)

    def reactionFrom(self, *molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """

        return self.weaver.reactionFrom(self, molecules)

    def reactionTo(self, *molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """

        return self.weaver.reactionTo(self, molecules)

    def setWeaver(self, weaver):
        """Internal - Should not be used outside of the framework.
        Sets the Weaver Object of this Circuit
        
        | *Args:*
        |     weaver: A weaver object that will be used
            """
        self.weaver = weaver

    def addAdvice(self, pointcut, adviceMethod, precedence=Advice.MINPRECEDENCE):
        """Declare a new advice in the aspect
        
        | *Args:*
        |     pointcut - The pointcut for the advice
        |     adviceMethod - The method to be executed at the pointcut
        |         should be a method bound to this aspect, with second parameter
        |         expecting a PointCutContext object
        |     precedence : integer (optional) - set precedence of advice,
        |         see notes on Advice.precedence
                
        | *Notes:*
        |     addAdvice constructs an Advice object. Further information thus
        |     can be found there.
        """

        self.adviceList.append(Advice(pointcut, adviceMethod, precedence))

    def addTypeAdvice(self, signature, typeaddition, name):
        """Declare a type advice in the aspect
        
        | *Args:*
        |     signature - The signature for the type advice
        |     typeaddition - The attribute / method to be added to the type
        |     name: str - The name for the typeaddition in the new type
            
        | *Notes:*
        |     addTypeAdvice constructs an TypeAdvice object.
        |     Further information thus can be found there.
        """

        self.typeAdviceList.append(TypeAdvice(signature, typeaddition, name, self))

    def addWeaverOutput(self, outputmethod):
        """Declare a new weaver output target 
        
        | *Args:*
        |     outputMethod - The method to be added to the WeaverOutput
            
        | *Raises:*
        |     InvalidWeaverOutputMethodError -
        |         If outputmethod is not a method or has wrong number of parameters
        |         (needs to accept self and a WeaverOutput reference)
        """

        if not inspect.ismethod(outputmethod) or len(inspect.getargspec(outputmethod)[0]) != 2:
            raise InvalidWeaverOutputMethodError("outputmethod must be a method with 2 parameters (self, WeaverOutput)")
        self.weaverOutputList.append(outputmethod)

    def getAdviceList(self):
        """Returns adviceList"""

        return self.adviceList

    def getTypeAdviceList(self):
        """Returns typeAdviceList"""

        return self.typeAdviceList

    def getWeaverOutputList(self):
        """Returns weaverOutputList"""

        return self.weaverOutputList

    @abstractmethod
    def mainAspect(self):
        """Entry point for an aspect, analogous to "main" in a program
        | Needs to be implemented by any sub class.
        | mainAspect will be called by the AOSB Weaver."""
        pass


class Weaver(object):
    """The "compiler" that weaves core concerns (circuits) and cross-cutting
    concerns (aspects) and creates a woven execution flow of parts
    
    | *Attributes:*
    |     partList - The current list of parts
    |     moleculeList - The current list of molecules
    |     beforeAndReplaceAdviceList - List of all before and replace advice to be woven
    |     afterAdviceList - List of all after advice to be woven
    |     partTypeAdviceList - List of all part type advice
    |     moleculeTypeAdviceList - List of all molecule type advice
    |     circuit - The main circuit
    |     aspects - The list of all aspects to be woven
    |     weaverOutput : WeaverOutput - The "compiled" result
    """


    class WeaverOutput(object):
        """Container for the woven result of the Weaver
        
        | *Attributes:*
        |     circuitName - the name of the circuit that was woven
        |     partList - finished ordered list of parts in the design
        |     moleculeList - list of all molecules in the design
        |     subcircuitList - list of all subcircuits (itself weaver outputs) (??)
        |
        """

        def __init__(self, circuitName, parentCircuit=None):
            self.circuitName = circuitName
            self.partList = []
            self.moleculeList = []
            self.wovenCircuitList = []
            self.parentCircuit = parentCircuit

        def __str__(self):
            result =  self.circuitName + ' Circuit Parts List:\n' + '+'.join(
                str(self.partList[i]) for i in range(len(self.partList)))

            result = result + "\nCircuit Molecule List: \n"
            result = result + self.molecules()

            result = result + "\nCircuit Molecule Reactions: \n"
            result = result + self.moleculeGraph()

            if len(self.wovenCircuitList) > 0:
                result = result + '\nSubCircuits:\n'
                for i in range(len(self.wovenCircuitList)):
                    result = result + "Sub Circuit "+str(i+1)+":\n"
                    result = result + str(self.wovenCircuitList[i])
                    result = result + '\n-----------\n'


            return result

        def molecules(self):
            return ', '.join(str(self.moleculeList[i]) for i in range(len(self.moleculeList)))

        def moleculeGraph(self):
            result = ''
            for i in range(len(self.moleculeList)):
                if(isinstance(self.moleculeList[i],Part) and self.moleculeList[i].scope != self):
                    continue
                for j in range(len(self.moleculeList[i].before)):
                    if (isinstance(self.moleculeList[i].before[j],Molecule)) or isinstance(
                            self.moleculeList[i].before[j], Part):
                        if(isinstance(self.moleculeList[i].before[j], Part)
                           and (self.moleculeList[i].before[j].scope == self)):
                            result += str(self.moleculeList[i].before[j].scope.circuitName) + '.' + str(self.moleculeList[i].before[j])\
                                      + '->' +\
                                      str(self.moleculeList[i].scope.circuitName) + '.'+ str(self.moleculeList[i])+'; '
                    else:
                        result += ','.join(str(self.moleculeList[i].before[j][k].scope.circuitName)+ '.' + str(self.moleculeList[i].before[j][k])
                                           for k in range(len(self.moleculeList[i].before[j])))\
                                  + '->' \
                                  + str(self.moleculeList[i].scope.circuitName)+ '.' +  str(self.moleculeList[i])+'; '

                for j in range(len(self.moleculeList[i].after)):
                    if((isinstance(self.moleculeList[i].after[j],Part) and self.moleculeList[i].after[j].scope != self)):
                        pass
                    else:
                        result += str(self.moleculeList[i].scope.circuitName) + '.' + str(self.moleculeList[i])\
                                  + '->' +\
                                  str(self.moleculeList[i].after[j].scope.circuitName) + '.' + str(self.moleculeList[i].after[j])+'; '

            return result

    def __init__(self, circuit, *aspects):
        """Sets of the weaver and compiles the design
        
        | *Attributes:*
        |     circuit : Circuit - The main circuit to be set
        |     *aspects : Aspect - list of aspects to be set
        
        | *Raises:*
        |     CircuitValueError - If circuit is invalid
        |     AspectValueError - If any aspect is invalid"""

        self.beforeAndReplaceAdviceList = []
        self.afterAdviceList = []
        self.partTypeAdviceList = []
        self.moleculeTypeAdviceList = []

        self.circuit = None
        self.aspects = []
        self.withinStack = []

        self.weaverOutput = None

        if not issubclass(circuit, Circuit):
            raise CircuitValueError("circuit must be a class of type Circuit.")

        self.circuit = circuit()
        self.circuit.setWeaver(self)

        for aspect in aspects:
            if not issubclass(aspect, Aspect):
                raise AspectValueError("all aspects must be classes of type Aspect.")
            self.aspects.append(aspect())

        self.weaverOutput = self.WeaverOutput(self.circuit.__class__.__name__)

        self.readAspectsConstructAdviceLists()

        self.sortAdviceList()

        # weaving is kicked off here
        self.currentWeaverOutput = self.weaverOutput
        self.circuit.mainCircuit()

        self.constructMoleculeListAndAddTypeAdvice()


    def sortAdviceList(self):
        """Internal - sorts the advice lists by precedence"""

        def sortKey(advice):
            return advice.precedence

        self.beforeAndReplaceAdviceList = sorted(self.beforeAndReplaceAdviceList, key=sortKey, reverse=True)
        self.afterAdviceList = sorted(self.afterAdviceList, key=sortKey)

    def constructMoleculeListAndAddTypeAdvice(self):
        pass

    '''    """Internal - constructs list of Molecules in design and adds their type advice"""

        # construct moleculeList - this should probably be done during weaving?
        for part in self.partList:
            if part.moleculeConnection:
                if not part.moleculeConnection in self.moleculeList:
                    
                    self.moleculeList.append(part.moleculeConnection)
                
                
                    
        # add type advice to molecules
        for molecule in self.moleculeList:
            for typeAdvice in self.moleculeTypeAdviceList:
                if typeAdvice.signature.match(molecule):
                    setattr(molecule,typeAdvice.name,typeAdvice.typeaddition)'''

    def readAspectsConstructAdviceLists(self):
        """Internal - initializes aspects and constructs all advice lists"""

        for aspect in self.aspects:
            aspect.setWeaver(self)
            aspect.mainAspect()
            for advice in aspect.getAdviceList():
                if advice.pointcut.operator == PointCut.AFTER:
                    self.afterAdviceList.append(advice)
                else:
                    self.beforeAndReplaceAdviceList.append(advice)

            for typeAdvice in aspect.getTypeAdviceList():
                if typeAdvice.isPartAdvice():
                    self.partTypeAdviceList.append(typeAdvice)
                else:
                    self.moleculeTypeAdviceList.append(typeAdvice)

            weaverOutput = self.weaverOutput  # to make weaverOutput visible to closure scope

            if len(aspect.getWeaverOutputList()) > 0:
                # this aspect wants to add a new weaver output target
                # TODO check if a similarly named output target already exists and raise 
                # useful exception (containing name of clashed method, name of advice)
                for outputMethod in aspect.getWeaverOutputList():
                    def outputMethodWrapper(method):
                        def outputMethodWrapperClosure(self):
                            return method(weaverOutput)

                        return outputMethodWrapperClosure

                    setattr(self.weaverOutput, outputMethod.__name__,
                            types.MethodType(outputMethodWrapper(outputMethod), aspect))

     # todo clean up, possible split lookup and creation
    def getMoleculeObject(self, scope, moleculeClass,createIt = False):
        for molecule in scope.moleculeList:
            if isinstance(molecule,moleculeClass):
                return molecule #TODO what about subclases?

        # molecule does not yet exist in the scope
        # (this means it's not been created or imported from another scope)
        if createIt == False:
            raise MoleculeValueError("Molecule " + moleculeClass.__name__ + " does not exist in this scope")
        else:
            result = moleculeClass()
            result.scope = scope
            result.additionStack = self.withinStack[:] #important: needs to be a deep copy!
            scope.moleculeList.append(result)
            self.addElemTypeAdvice(result)
            return result

    def addElemTypeAdvice(self, elem):
        """Internal - adds type advice to a part"""

        list = self.partTypeAdviceList
        if isinstance(elem,Molecule):
            list = self.moleculeTypeAdviceList

        for typeAdvice in list:
            if typeAdvice.signature.match(elem):
                # TODO check for name clashes...
                if inspect.ismethod(typeAdvice.typeaddition):
                    def adviceMethodWrapper(method):
                        def adviceMethodWrapperClosure(self):
                            return method(elem)

                        return adviceMethodWrapperClosure

                    setattr(elem, typeAdvice.name,
                            types.MethodType(adviceMethodWrapper(typeAdvice.typeaddition), typeAdvice.aspect))
                else:
                    setattr(elem, typeAdvice.name, typeAdvice.typeaddition)

    def runBeforeAndReplaceAdvice(self, callingObject, part):
        """Internal - runs Before and Replace Advice of a part
        
        | *Args:*
        |     part - The part for which the advice matches
        |     callingObject - the circuit or advice which added the part
        
        | *Returns:*
        |     "continue"
        |     True - If part should still be added
        |     False - If part has been replaced and remaining after advice executed
        """

        for advice in self.beforeAndReplaceAdviceList:
            if len(self.currentWeaverOutput.partList) > 0:
                part.setBeforePart(self.currentWeaverOutput.partList[-1])

            if advice.pointcut.match(part):
                if (advice.pointcut.operator == PointCut.REPLACE):
                    numberOfMatchingParts = advice.pointcut.signature.numberOfMatchingParts(part)
                    while numberOfMatchingParts > 1:
                        self.currentWeaverOutput.partList.pop()
                        numberOfMatchingParts -= 1
                    advice.adviceMethod(PointCutContext(self.withinStack, part))

                    self.runAfterAdvice(callingObject, part, advice.precedence)
                    return False
                advice.adviceMethod(PointCutContext(self.withinStack, part))

        return True

    def runAfterAdvice(self, callingObject, part, precedence=Advice.MINPRECEDENCE):
        """Internal - runs After Advice of a part
        
        | *Args:*
        |     part - The part for which the advice matches
        |     callingObject - the circuit or advice which added the part
        |     precedence - Only run advice with a precedence greater or equal this
        |         Used if replacement advice has been executed
        
        """

        # after advice are ordered by ascending precedence
        #the one with the lowest priority is executed
        #if precedenceLevel is set by a replacement advice, we will only
        #execute those advice with >= precedence
        for advice in self.afterAdviceList:
            if (advice.precedence >= precedence) and (advice.pointcut.match(part)):
                advice.adviceMethod(PointCutContext(self.withinStack, part))

    def addReaction(self, callingObject, fromMoleculeList, toMoleculeList):
        realFromMoleculeList = []
        realToMoleculeList = []
        for element in fromMoleculeList:
            realFromMoleculeList.append(self.getMoleculeObject(self.currentWeaverOutput,element))
            # todo needs to be adapted to actually search scope ??? or report error
        for element in toMoleculeList:
            realToMoleculeList.append(self.getMoleculeObject(self.currentWeaverOutput,element,True))

        # todo each element in fromMoleculeList does not need toMoleculeList in its before...
        # need to show that here is a gate...

        # if both are just one molecule, can hook them up
        if (len(realFromMoleculeList) == 1) and (len(realToMoleculeList) == 1):
            realFromMoleculeList[0].after.append(realToMoleculeList[0])
            realToMoleculeList[0].before.append(realFromMoleculeList[0])
        elif (len(realFromMoleculeList) > 1) and (len(realToMoleculeList) == 1):
            realToMoleculeList[0].before.append(realFromMoleculeList)
        elif (len(realFromMoleculeList) == 1) and (len(realToMoleculeList) > 1):
            realFromMoleculeList[0].after.append(realToMoleculeList)



    class MoleculeReactionTo:
        def __init__(self,molecules):
            self.molecules = molecules

    class MoleculeReactionFrom:
        def __init__(self,weaver,callingObject, molecules):
            self.weaver = weaver
            self.callingObject = callingObject
            self.fromMolecules = molecules

        def __rshift__(self,other):
            # todo check type of other, should be MoleculeReactionTo
            self.weaver.addReaction(self.callingObject,self.fromMolecules,other.molecules)


    def reactionFrom(self, callingObject, molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """
        return self.MoleculeReactionFrom(self,callingObject,molecules)


    def reactionTo(self, callingObject,  molecules):
        """Used to add parts to the circuit by passing them on to the AOSB Weaver

        | *Args:*
        |     part: The part to be added to the circuit

        """
        return self.MoleculeReactionTo(molecules)

    def importMolecule(self, callingObject, molecule):
        # means that we are getting a molecule from the next larger scope
        # todo check that we are not in outest scope
        # todo rethink - why should getMOlecule... take an object?? it's creating them afterall... should take class name?
        # todo no? maybe?
        importedMoleculeParent = self.getMoleculeObject(self.currentWeaverOutput.parentCircuit,molecule)

        # case 1: molecule already exists in scope and needs to be merged
        for mol in self.currentWeaverOutput.moleculeList:
            if isinstance(mol,molecule):
                mol.before.append(importedMoleculeParent)
                importedMoleculeParent.after.append(mol)
                return

        #otherwise:
        #case 2: molecule doesn't exist yet
        self.createMolecule(callingObject,molecule)
        # now we can add it via case 1:
        self.importMolecule(callingObject,molecule)


    def exportMolecule(self, callingObject, molecule):
        # todo error checking
        # e.g. if molecule already exists (but possible differnet instances?? warn of clash?
        # if parent doesn't exist
        exportMolecule = self.getMoleculeObject(self.currentWeaverOutput,molecule)

        # case 1: the molecule already exists in the outer scope and needs to be added to its graph
        for mol in self.currentWeaverOutput.parentCircuit.moleculeList:
            if isinstance(mol,molecule):
                mol.before.append(exportMolecule)
                exportMolecule.after.append(mol)
                return
        #otherwise:
        #case 2: molecule doens't exist yet in outer scope
        self.getMoleculeObject(self.currentWeaverOutput.parentCircuit,molecule,True)
        self.exportMolecule(callingObject,molecule)


    def createMolecule(self, callingObject, molecule,scope = None):
        self.getMoleculeObject(self.currentWeaverOutput,molecule,True)
        # todo, if create is set to true, should actually report an error if molecule already exists.


    def addCircuit(self, callingObject, circuit):
        # todo check circuit
        circuitObject = circuit()
        newCircuit = self.WeaverOutput(circuitObject.__class__.__name__,self.currentWeaverOutput)
        self.currentWeaverOutput.wovenCircuitList.append(newCircuit)
        self.currentWeaverOutput = newCircuit

        # new sub circuit is set up. no we'll weave it

        circuitObject.setWeaver(self)
        circuitObject.mainCircuit()

        self.currentWeaverOutput = self.currentWeaverOutput.parentCircuit

    def addPart(self, callingObject, part):
        """Called by circuit or aspect to add a part in the execution flow
        
        | *Args:*
        |     callingObject - the circuit or aspect calling this
        |     part - The part supposed to be added
        
        """

        part = checkIfTypeReturnInstance(part)
        # part.namespace = self.currentWeaverOutput
        part.scope = self.currentWeaverOutput


        self.withinStack.append(callingObject)
        part.additionStack = self.withinStack[:]
        continueAddingPart = self.runBeforeAndReplaceAdvice(callingObject, part)


        # before adding the part, add any type advice



        if continueAddingPart:
            part.weave(self)
            # TODO clean this up, the namespace issue...


            self.addElemTypeAdvice(part)


            # If the part is a composite, we unpack it here. otherwise we finally add the part  
            if isinstance(part, Circuit) == False:

                #add before / after information
                if len(self.currentWeaverOutput.partList) > 0:
                    self.currentWeaverOutput.partList[-1].setAfterPart(part)
                    part.setBeforePart(self.currentWeaverOutput.partList[-1])

                self.currentWeaverOutput.partList.append(part)

            else:
                # the circuit is not initialized with a weaver
                part.setWeaver(self)
                part.mainCircuit()

            self.runAfterAdvice(callingObject, part)

        self.withinStack.pop()


    def output(self):
        """Returns the weaverOutput Element"""
        return self.weaverOutput


class MoleculeValueError(ValueError):
    """Molecule was expected, but something else was given"""
    pass


class PartInitializationError(Exception):
    """Part was expected, but something else was given"""


class PartValueError(ValueError):
    """Part was expected, but something else was given"""


class InvalidSymbolNameError(ValueError):
    """A symbol name (string) was not correctly formatted"""
    pass


class SymbolExistsWarning(UserWarning):
    """A Part or Molecule declaration uses an existing name"""
    pass


class InvalidSignatureError(ValueError):
    """A signature (PointCut / Part / Molecule) is incorrectly formatted or typed"""
    pass


class InvalidPointCutExpressionError(ValueError):
    """There is an error in a Point Cut expression"""
    pass


class InvalidPointCutOperatorError(ValueError):
    """An unknown or illegal operator was used for the point cut"""
    pass


class InvalidPointCutError(ValueError):
    """object of type PointCut expected, but something else given"""
    pass


class InvalidAdviceMethodError(ValueError):
    """method with 2 parameters expected, but something else given"""
    pass


class InvalidWeaverOutputMethodError(ValueError):
    """ with 2 parameters expected, but something else given"""
    pass


class PrecedenceOutOfRangeError(ValueError):
    """An unknown or illegal operator was used for the point cut"""
    pass


class AspectValueError(ValueError):
    """An Aspect class was expected, but something else given"""


class CircuitValueError(ValueError):
    """A Circuit class was expected, but something else given"""