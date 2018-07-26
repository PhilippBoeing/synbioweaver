# location of synbioweaver
export PYTHONPATH=$PYTHONPATH:../../

# location of cuda-sim
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/cuda-sim
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/abc-sysbio

python run_constGFP_inference.py 
