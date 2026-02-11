## Anaconda package for private channel

```
# Build environment
conda create -n conda-build -c conda-forge python=3.12 conda-build anaconda-client boa -y
conda activate conda-build

# Version from project
export MUTME_VERSION="$(python -c 'import tomllib; print(tomllib.load(open("pyproject.toml","rb"))["project"]["version"])')"
echo "$MUTME_VERSION"

# From repo root
conda-build -c conda-forge -c bioconda distr/conda/recipe

# Find build artifact
conda-build --output distr/conda/recipe

# Upload to channel after login
anaconda upload --user <user_channel> "$(conda-build --output distr/conda/recipe)"
```

## PyPI package

```
# Install development and test dependencies
pip install .[dev,test]

# From repo root to run tests
pytest

# Build package
python -m build

# Upload to PyPI after token is configured
python -m twine upload dist/mutme-*
```