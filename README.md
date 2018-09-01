synbioweaver
=============

### Installation
To use the package, the python module `synbioweaver` must be on the PYTHONPATH.

### Dependencies
The usual python dependencies are required for running the base functionality: `numpy`, `scipy`, `matplotlib`. If you would like to export designs into sbml then you also need `libsbml`. The more involved examples require additional tools and are outlined below

### Documentation
An overview and tutorial is provided at the following location
<http://synbioweaver.readthedocs.io/en/latest/>

### Examples
A number of simple examples to illustrate how the package works are included in the folder examples/tutorial. More involved examples demonstrating the flexibility and applicability of the framework are given in the other folders in examples. These include:

**Weaving together gene circuits**  
examples/weaving-circuits  

**Generation of transfer curves based on part libraries**  
examples/transfer-curves  

**Specifying different gene circuit designs and generating models**  
examples/model-generation  
The circuits include: 
* simple implementations of constitutive and inducible GFP expression
* the lux-AHL system
* the sender-receiver system found in Basu et. al. Spatiotemporal control of gene expression with pulse-generating networks (2004) <http://www.pnas.org/content/101/17/6355>
* the toggle-switch system found in Litcofsky et. al. Iterative plug-and-play methodology for constructing and modifying synthetic gene networks (2012) <https://www.nature.com/articles/nmeth.2205>
* the oscillator from Stricker et. al. A fast, robust and tunable
  synthetic gene oscillator (2008) <https://www.nature.com/articles/nature07389>
* the system described in Ceroni et. al. (2018)
  <https://www.nature.com/articles/nmeth.4635> which uses combination
  of guide RNA and dCas9 to provide burden driven negative feedback
* a generic design for a whole-cell biosensor that makes use of an
  expressed two component system to drive expression of a reporter

**Multicellular logic gates**  
examples/logic-gates  
This example is based on the XOR gate described in Tamsir
et. al. Robust multicellular computing using genetically encoded NOR
gates and chemical ‘wires’ (2011)
<https://www.nature.com/articles/nature09565>. It also shows how
aspects can be used to generate input-output transfer curves for different circuits
using characterisation data.

**Rule-based modelling**  
examples/rule-based-model   
A design for the repressilator is converted into a rule-based model specified in the Kappa language and then subsequently simulated stochastically using `KaSim`. This example requires `KaSim` version 3.5 to be installed in /usr/local/share. This can be downloaded as a binary from <https://github.com/Kappa-Dev/KaSim/releases>

**Automated model generation and simulation using cuda-sim**  
examples/context-simulation  
For GPU based biochemical network simulation requires a suitable Nvidia GPU device, `CUDA`, `PyCUDA` must be installed and `cuda-sim` must be on PYTHONPATH

**Bayesian inference using ABC-SysBio**  
examples/characterisation  
In this example the OD and fluorescence of a simple constitutive GFP circuit are simulated over time. The resultant data are fit to an ODE model of the system using Bayesian statistics. This currently uses `cuda-sim`. In addition it requires the module `abcsysbio` to be on the PYTHONPATH. 

**Reverse engineering context dependence**  
examples/context-inference  
This example is inspired by Catanach et. al. Context Dependence of Biological Circuits (2018) <https://www.biorxiv.org/content/early/2018/07/03/360040>. Multiple models are generated of a circuit under two contexts. Bayesian model selection is applied to formally test whether the systems have the same parameters. 

