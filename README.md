# Polycat: Framed Polytopes And Higher Categories

This repository accompanies the paper **"Framed Polytopes and Higher Categories"** by Guillaume Laplante-Anfossi, Anibal M. Medina-Mardones, and Arnau Padrol. It contains SageMath code for computing **cellular strings** and **loops** on framed polytopes.

## Contents

- `polycat.sage`: SageMath script containing core functions for computations related to framed polytopes.
- `polycat.ipynb`: Jupyter notebook with examples illustrating how to use the code in the context of the paper.

## Requirements

- [SageMath](https://www.sagemath.org/) (version ≥ 9 recommended)

## Usage

To explore examples and run the computations, open the notebook:

```bash
sage -n jupyter
```

Then open `polycat.ipynb` in your browser.

The notebook uses:

```python
load("polycat.sage")
```

So make sure that `polycat.sage` is in the same directory as `polycat.ipynb`, or adjust the path in the `load()` command accordingly.


## License

The code is licensed under the [GNU General Public License v3 or later](https://www.gnu.org/licenses/gpl-3.0.html).
