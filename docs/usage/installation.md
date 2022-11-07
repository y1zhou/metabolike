# Installation

!!! note
This package is still under very active development, and **anything** may change at **any moment**. [Semantic versioning](https://semver.org/) will be followed starting from release v1.0.0.

## Releases from PyPI

To install the latest stable version of the package, run the following command:

```bash
pip install metabolike
```

If you are being a bit more adventurous, the pre-release version is also available:

```bash
pip install --pre metabolike
```

## Latest version from GitHub

As the package is still in its early stages, it is recommended to get the latest version from GitHub:

```bash
git clone --depth 1 https://github.com/y1zhou/metabolike
```

To install the package, run the following command:

```bash
pip install . # (1)
```

1. Note that you should `cd` into the directory first, i.e. `#!bash cd metabolike`.

### Contributing

If you would like to contribute to `metabolike`, some extra dependencies including [pytest](https://docs.pytest.org/), [MkDocs](https://www.mkdocs.org/), [black](https://github.com/psf/black), and [pre-commit](https://pre-commit.com/) are needed:

```bash
pip install -e .[dev]
```
