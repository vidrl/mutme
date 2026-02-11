## Pre-commit hooks

See `.pre-commit-config.yaml` for available hooks.

```
# Install package with dev-dependencies
pip install .[dev]

# Install pre-commit hook for Ruff
pre-commit install

# Custom hook script `bump_version_commit_msg.sh`
# Install version bump on feat/fix/hotfix semantic version commit messages
pre-commit install --hook-type commit-msg
```