import remote_runner
from remote_runner import Pool, LocalSlurmWorker, Task
from pathlib import Path
from glob import glob

workers = [
    # worker which sync periodically and at the end of task
    LocalSlurmWorker(
        remote_user_rc = """
unset PYTHONPATH
source ~/cov2/venv/bin/activate
source /opt/amber20/amber.sh
""",
        sbatch_args = [
            #"--gpus=1", # Select 1 arbitrary gpu
            "--gpus=rtx_3080:1", # Alternative: select 1 rtx_2080_ti card
            "--cpus-per-task=2", # cpu here means threads 
            "--ntasks=1", # task is equivalent to number of processes 
                          # ntasks >=1 makes sense only if MPI is involved
            #"--nodelis=bionmr-mom-006"
        ])
    for i in range(3) # number of parallel runners / jobs
]

mutants = ["wt+mp1", "alpha+mp1", "delta+mp1", "delta_plus+mp1", "omicron+mp1"]

tasks = [
    Task.load(Path(f"../sample/{mutant}/state.dill") )
    for mutant in mutants
]
# print(tasks)

#states = sorted(glob("lcb1_meta12/state.dill"))
#print(states)
#tasks = [Task.load(Path(state)) for state in states]
Pool(workers).run(tasks)

