Before everything starts, please make sure Conda has been properly installed.
The command below will help to install all the software included in the core environment:

`conda env create --file=sgp_env.yml`

Then try activating it:

```
source activate sgp_env
```

Although `source activate` is functionally interchangeable with `conda activate`, we still recommend to test with the former. It will be required by snakemake, so these commands could also help to ensure. 

Then we can try to ensure the environment does work:

```
which snakemake
which python
```

If these commands show some paths one under your Conda environment folder (e.g. ~/miniconda3/envs/sgp_env/bin/snakemake and ~/miniconda3/envs/sgp_env/bin/python, respectively), the environment already be successful. 

When they look fine, we can turn the environment off: 

```
source deactivate
```

Please also note that the environment should not have changed any original software settings, so everything, such as environment variables or software versions, should be same as what they were before the installation was done. The software included in the environment is available only when the environment has been activated. 

