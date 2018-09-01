# location of synbioweaver
export PYTHONPATH=$PYTHONPATH:../../

# location of cuda-sim
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/cuda-sim
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/abc-sysbio


rm -r res_abc; python -u run_indGFP_context_inference.py
