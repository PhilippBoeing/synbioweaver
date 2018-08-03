# location of synbioweaver
export PYTHONPATH=$PYTHONPATH:../../

# location of cuda-sim and abc-sysbio
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/cuda-sim
export PYTHONPATH=$PYTHONPATH:/home/cbarnes/soft/abc-sysbio

python run_toggleOscillator_combined.py > log.txt
python run_toggleOscillator_coupA_AraC.py >> log.txt
python run_toggleOscillator_coupA_GFP.py >> log.txt
