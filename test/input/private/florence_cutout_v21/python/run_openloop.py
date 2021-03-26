# Purpose is to advance openloop to get some restarts.
# Restart on 2018-08-01 and run until 2018-08-15

import datetime
import deepdiff
import os
import pathlib
import pickle
import time
import wrfhydropy

florence_path = pathlib.Path('/glade/scratch/jamesmcc/tmp_florence_cutout_run_openloop')

if not florence_path.exists():
    os.mkdir(florence_path)

#config = 'dart_nwm_channel-bucket-only'
config = 'nwm_ana'
domain_path = pathlib.Path('/glade/work/jamesmcc/domains/private/florence_cutout_v21/')
hydro_namelist_config_file= domain_path / 'hydro_namelists.json'
hrldas_namelist_config_file= domain_path / 'hrldas_namelists.json'
compile_options_config_file= domain_path / 'compile_options.json'
    
comp_path = florence_path / 'compile'
# for a fresh compile, remove comp_path from disk.
if not comp_path.exists():

    model = wrfhydropy.Model(
        '/glade/u/home/jamesmcc/WRF_Hydro/wrf_hydro_nwm_public_2/trunk/NDHMS',
        model_config=config,
        compiler='ifort',
        hydro_namelist_config_file=hydro_namelist_config_file,
        hrldas_namelist_config_file=hrldas_namelist_config_file,
        compile_options_config_file=compile_options_config_file
    )

    ## make sure your modules match the compiler first...
    model.compile(comp_path)

else:
    model = pickle.load((comp_path / 'WrfHydroModel.pkl').open('rb'))
    model.model_config = config  # cheat/hack

domain = wrfhydropy.Domain(
    domain_top_dir=domain_path,
    domain_config=config
)


domain_restart = domain.hrldas_namelist_patches['noahlsm_offline']['restart_filename_requested']
model_start_time = datetime.datetime(2018, 8, 1)

model_end_time = model_start_time + datetime.timedelta(hours=24*75)
job = wrfhydropy.Job(
    exe_cmd='mpirun -np 2 ./wrf_hydro.exe',
    job_id='flo_cut',
    restart=True,
    model_start_time= model_start_time,
    model_end_time=model_end_time
)

# scheduler = wrfhydropy.schedulers.PBSCheyenne(
#     account='NRAL0017',
#     email_who='jamesmcc@ucar.edu',
#     email_when='abe',
#     nproc=36,
#     nnodes=1,
#     ppn=36,
#     queue='premium',
#     walltime='00:10:00'
# )

sim = wrfhydropy.Simulation()
sim.add(model)
sim.add(domain)
sim.add(job)
# sim.add(scheduler)
# # -------------------------------------------------------
# # Simulation Object
florence_run = florence_path / 'run_openloop'
os.mkdir(florence_run)
os.chdir(florence_run)
sim.compose()
sim.run()

