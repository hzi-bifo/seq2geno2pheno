Before everything starts, please make sure Conda has been properly installed.
The command below will help to install all the software included in the core environment:

`conda env create --file=sgp_env.yml`

Then try activating it:

```
source activate sgp_env
which snakemake
which python
source deactivate
```

If the command above shows some paths one under your Conda environment folder (e.g. ~/miniconda3/envs/sgp_env/bin/snakemake and ~/miniconda3/envs/sgp_env/bin/python, respectively), the environment already be successful. 
Please also note that the environment should not have changed any original software settings such as python version or available packages. The software included in the environment is available only when the environment has been activated. 
