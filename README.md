# z-quantum-solid-state

[![codecov](https://codecov.io/gh/zapatacomputing/z-quantum-solid-state/branch/master/graph/badge.svg?token=UPHNS662NE)](https://codecov.io/gh/zapatacomputing/z-quantum-solid-state) 

## What is it?

`z-quantum-solid-state` is a module for solid state calculations to be used with [Orquestra](https://www.zapatacomputing.com/orquestra/) – platform for performing computations on quantum computers developed by [Zapata Computing](https://www.zapatacomputing.com).


## Usage

### Workflow
In order to use `z-quantum-solid-state` in your workflow, you need to add it as a resource:

```
resources:
- name: z-quantum-solid-state
  type: git
  parameters:
    url: "git@github.com:zapatacomputing/z-quantum-solid-state.git"
    branch: "master"
```

and then import in a specific step:

```
- - name: my-task
    template: template-1
    arguments:
      parameters:
      - param_1: 1
      - resources: [z-quantum-solid-state]
```

Once it's done you can:
- use any template from `templates/` directory
- use tasks which import `zquantum.solidstate` in the python code (see below).

### Python

Here's an example how to do use methods from `z-quantum-solid-state` in a python task:

```python
from zquantum.solid_state.fermi_hubbard import calculate_exact_density_of_energy_for_2_D_fermi_hubbard
energy_density = calculate_exact_density_of_energy_for_2_D_fermi_hubbard(4,2,2,2)
```

Even though it's intended to be used with Orquestra, you can also use it as a standalone python module.
In order to install it run `pip install .` from the `src` directory.


## Development and contribution

You can find the development guidelines in the [`z-quantum-core` repository](https://github.com/zapatacomputing/z-quantum-core).

### Running tests

In order to run tests please run `pytest .` from the main directory.
