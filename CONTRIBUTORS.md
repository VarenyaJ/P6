# Contributors

Instructions for developers who want to contribute to the project.

We recommend to check out the repository.

```shell
git clone https://github.com/VarenyaJ/P6.git
cd P6
```

Then, setup the environment. You can use `venv` or `conda`.

Here is an example for conda:

```shell
conda env create -f requirements/environment.yml -y
```

or use PIP directly:

```shell
python3 -m pip install -r requirements/requirements.txt -r requirements/requirements_test.txt .
```

This will install `P6` along with the dependencies needed for the development.

Verify that the tests pass by running:

```
pytest
```
