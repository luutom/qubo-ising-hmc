# qubo_ising_hmc

## Description

## Install
Install via pip
```bash
pip install [--user] -r requirements.txt
pip install [-e] [--user] .
```
If you want to specify a specific CC file, you can install as (for example, a homebrewed gcc)
```
CC=/usr/local/Cellar/gcc/10.2.0/bin/c++-10  pip install -e .
```


## Run

See the `notebooks/example-usage.ipynb`.

To run the notebook, you also have to pip install [`quantum_linear_programming`](https://github.com/cchang5/quantum_linear_programming)
```bash
cd quantum_linear_programming
pip install [-e] [--user] .
```

## License
See [LICENSE](LICENSE)
