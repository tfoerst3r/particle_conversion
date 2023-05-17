---
tags: #1Dmodel, #coal-conversion, #python
title: One Dimension Particle Conversion Model
---

<h1 align="center">One Dimension Particle Conversion Model</h1>
<p align=center> 
This model simulates the conversion process, based on porous carbonaceous materials, like coal.
</p>

# About the Project (optional)

The project simulates different model approaches for the conversion process of porous carbonaceous particle. The conversion process depends on a variety of property and conditions, like overall pressure.

<!------------------->
<!----- Chapter ----->
<!------------------->
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








> (If applicable) add step-by-step instructions on how to install and run the project 

## Usage

> How to use your project? 
> How to handle everyday use cases? --> this can also be done on a separate documentation page






## Roadmap (optional)

> Give a glimpse of features that are on the roadmap

## Contributing (optional)

> Put other information here, like where to get help and find more project-related information, 
> e.g., community forums, discord channels, Trello boards 

We are accepting merge requests, bug reports, and feature requests. 
Please see the [CONTRIBUTING.md](CONTRIBUTING.md) for the complete guide. 

## License

Distributed under the MIT License. See `LICENSE.md` for more information.

-----

> There are many more sections that can be included here, like 
>
> * Contributors 
> * Contact information 
> * FAQ


# Tray

Possible adaption:

- self.dXdt = 0.5 * self.dXdt
    - bader pS approach
    - fluent pS approach

- self.Deff = 0.5 * self.Deff


