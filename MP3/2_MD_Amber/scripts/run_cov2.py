import remote_runner
from remote_runner import Pool, LocalSlurmWorker, Task
from pathlib import Path
from glob import glob

workers = [
    # worker which sync periodically and at the end of task
    LocalSlurmWorker(
        remote_user_rc = """
unset PYTHONPATH
source ~/cov2/venv_cov2/bin/activate
source /opt/amber20/amber.sh
""",
        sbatch_args = [
            "--gpus=rtx_2080_ti",           # Select gpu type from your machine
            "--cpus-per-task=2",            # cpu here means threads 
            "--ntasks=1",                   # task is equivalent to number of processes 
                                            # ntasks >=1 makes sense only if MPI is involved
            "--nodelis=bionmr-mom-003"      # scpecify machine name
        ])
    for i in range(7) # use to set a number of parallel runners / jobs
]

# mutants for simulation
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

tasks = [
    Task.load(Path(f"./{mutant}/state.dill") )
    for mutant in mutants
]

Pool(workers).run(tasks)

