# pmvalsampsize

Sample size for external validation of a prediction model

## Installation

```bash
$ pip install pmvalsampsize
```

## Usage

Computes the minimum sample size required for the external validation of an existing multivariable prediction model using the criteria proposed by [Archer (2020)](doi:10.1002/sim.8766) and [Riley (2021)](doi:10.1002/sim.9025)

```python
from pmvalsampsize.pmvalsampsize import *

samplesize = pmvalsampsize(type="b", cstatistic=0.8, graph=True, lpbeta=(1.33, 1.75), prevalence=0.43, noprint=True) # change options to meet youor requirements
plt.show()
summary(samplesize)
```

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`pmvalsampsize` was created by Rebecca Whittle. It is licensed under the terms of the GNU General Public License v3.0 license.

## Credits

`pmvalsampsize` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
