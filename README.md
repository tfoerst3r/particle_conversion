<!--
SPDX-FileCopyrightText: 2023 Thomas Förster

SPDX-License-Identifier: CC-BY-4.0
-->

<h1 align="center">One Dimension Particle Conversion Model</h1>
<p align=center> 
This model simulates the conversion process, based on porous carbonaceous materials, like coal.
</p>


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# About the Project

The project simulates different model approaches for the conversion process of porous carbonaceous particle. The conversion process depends on a variety of property and conditions, like overall pressure.


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Getting Started

## Prerequisites

- To work with the package you need to have at least Python 3.10 installed.

- This project uses the Python dependency manager `poetry` [(Installation guide for poetry)][poetry-install]

    Example command line for Linux/Ubuntu:

    ~~~
    $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3 -
    ~~~
    
    or

    ```
    $ pip install poetry
    ```

## Installation

### Using the environment distributed with this model

**Step 1. Clone the repository.**

- via ssh (first, create [SSH key][ssh-key])

    ~~~
    $ git clone --depth=1 --branch=main git@github.com:tfoerst3r/particle_conversion.git 
    ~~~ 

- or via https

    ~~~
    $ git clone --depth=1 --branch=main https://github.com/tfoerst3r/particle_conversion.git
    ~~~ 

Alternatively, you can download source code directly via the download option on top of the repository page.


**Step 2. Install required dependencies.**

To build the package run the following in the root folder of the project:

```
$ poetry install
```

### Using any other python environments

Poetry can build, so called wheels, which can be installed in any other python environment.

**Step 1. Create the wheel or download the wheel**

You can download any wheel via any desired CI job called "build_wheel_package". Or build from scratch via the following step:

```
$ poetry build --format wheel
```

This will create a wheel file `particle-conversion-*.whl` in the folder `./dist`.

**Step 2. Installing the wheel**

- First you need to activate your desired python environment, e.g. `$ source /path/to/your/env/bin/activate`. And then install your the wheel.

```
$ pip install particle_conversion-*.whl
```

- With `$ which pconv` you can check if the executable is in the desired environment.


## Uninstall the module:

For uninstalling `pconv` in your current environment use:

~~~
$ pip uninstall pconv
~~~


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Usage


## Prerequisites

First you need to activate the environment you want to work in. Either use

- the provided poetry environment

    ```
    $ poetry serve
    ```

- or any other environment where you installed the provided wheel

    ```
    $ source /path/to/your/env/bin/activate
    ```

In both versions an executable `pconv` should be available. Alternatively you can run the CLI directly via:

```
$ poetry run pconv [OPTIONS]
```


## Usage of the `pconv` script

**Syntax**

```
$ pconv --config [CONFIG FILE] --output [OUTPUT FILE]
```

For more help type:

```
$ pconv --help
```

**Step 1. Choose the appropriate model settings and adapted the corresponding configuration file.**

- Choose the appropriate growth model (for more details see details in [the model file](./docs/conversion_models.md)) 
- Adapt the corresponding model configuration `*.toml` file (see base example in [config file](./config/).

**Step 2. Run the model with the desired output path.**


## Examples

**default model**

```
$ pconv --config ./config/base.toml --output ./output/base_output.csv
```


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Contributing

**We welcome your contribution!**

The repository is still under development and any feedback, suggestions, technical contributions are highly welcome.

General Options:

- open an issue, if you have improvement ideas
- fork the project and contribute via merge request against the main branch of this repository


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# License

Please see the file [LICENSE.md](./LICENSE.md) for further information about how the content is licensed.

<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Citation

If you are using this code in your publications, please refer to [DOI:10.](https://doi.org/) for citation, or cite as:

<small>
Thomas Förster. (2023). <i>Particle Conversion Solver: Testing tool for different conversion models via custom predefined boundary conditions.</i> Zenodo. <a href="https://doi.org/">https://doi.org/</a>
</small>


<!-- Literature -->

<!-- Links -->

