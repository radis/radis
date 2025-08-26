# GitHub Actions Workflows

This directory contains GitHub Actions workflows that replace the Travis CI configuration.
The workflows are designed to provide the same functionality as the original `.travis.yml` file.

## Workflows Overview

### 1. CI Workflow (`ci.yml`)
**Triggers:** Push/PR to `master` and `develop` branches

**Jobs:**
- **Test Matrix:** Tests on multiple Python versions (3.10, 3.11) and Ubuntu 20.04
- **Extended Testing:** Additional Python versions (3.8, 3.12) and macOS testing (master branch only)
- **Fast Tests:** Quick test suite with database exclusions
- **Long Tests:** Comprehensive test suite (skipped on macOS)
- **Coverage:** Uploads coverage to Codecov (Python 3.10 on Ubuntu only)

**Key Features:**
- Uses Micromamba for conda environment management
- Headless display support with xvfb on Linux
- Conditional job execution based on branch
- Caching for faster builds

### 2. Code Quality Workflow (`lint.yml`)
**Triggers:** Push/PR to `master` and `develop` branches

**Jobs:**
- **Lint:** Runs pre-commit hooks for code quality checks
- **Error Reporting:** Shows detailed logs on failure

### 3. Deployment Workflow (`deploy.yml`)
**Triggers:** Push to `master` branch, manual dispatch

**Jobs:**
- **Deploy:** Builds and publishes to PyPI
- **Post-Deploy Testing:** Installs from PyPI and runs tests to verify

**Security:**
- Uses GitHub environment for deployment protection
- Requires `PYPI_API_TOKEN` secret

## Required Secrets

Set these in your repository settings under Settings → Secrets and variables → Actions:

1. **`PYPI_API_TOKEN`**: PyPI API token for package publishing
2. **`CODECOV_TOKEN`**: Codecov token for coverage reporting (optional but recommended)

## Environment Setup

Create a GitHub environment named `deployment` in your repository settings for additional protection of the deployment workflow.

## Migration Notes

### Key Differences from Travis CI:

1. **Matrix Strategy:** GitHub Actions uses a more flexible matrix strategy
2. **Conditional Jobs:** Uses `if` conditions instead of Travis stages
3. **Shell Handling:** Explicit shell configuration for cross-platform compatibility
4. **Caching:** Built-in environment caching for faster builds
5. **Security:** Uses GitHub environments and secrets for secure deployment

### Equivalent Functionality:

| Travis CI Stage | GitHub Actions Workflow | Notes |
|----------------|-------------------------|-------|
| `test` | `ci.yml` main jobs | Standard testing on Ubuntu |
| `lint` | `lint.yml` | Code quality checks |
| `deploy` | `deploy.yml` | PyPI deployment |
| `Py38`, `Py312` | `ci.yml` matrix | Additional Python versions |
| `MacOS_py310` | `ci.yml` matrix | macOS testing |

### Test Markers Preserved:

All pytest markers from the original Travis configuration are maintained:
- `fast` vs `not fast`
- Database exclusions (`needs_db_*`)
- CUDA exclusions (`needs_cuda`)
- Large database exclusions (`download_large_databases`)

## Troubleshooting

### Common Issues:

1. **Micromamba Setup:** If environment creation fails, check the `environment.yml` syntax
2. **Test Failures:** Ensure all required system dependencies are installed
3. **Deployment Failures:** Verify PyPI token permissions and package version conflicts
4. **Coverage Upload:** Check Codecov token and repository configuration

### Debug Tips:

- Use `workflow_dispatch` trigger for manual testing
- Add `debug: true` to micromamba setup for verbose logging
- Check job logs for detailed error messages
- Verify secrets are properly configured

## Performance Optimizations

- **Caching:** Environment caching reduces setup time
- **Parallel Jobs:** Matrix jobs run in parallel
- **Conditional Execution:** Skip unnecessary jobs on feature branches
- **Artifact Management:** Automatic cleanup of build artifacts
