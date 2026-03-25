# Conda Packaging for MADNESS

Conda packages are built and uploaded to [anaconda.org/m-a-d-n-e-s-s](https://anaconda.org/m-a-d-n-e-s-s) automatically via GitHub Actions.

Users install MADNESS with:
```
conda install -c m-a-d-n-e-s-s madness
```

## Automated Workflow

The workflow is defined in `/.github/workflows/conda-deploy.yml`. It builds a minimal MADNESS (no MPI, no libxc) for Linux and macOS and uploads the conda package.

It runs automatically when:
- The CI tests pass on `master` (triggered by `workflow_run` after the "Linux/MacOS Build" workflow succeeds)
- A version tag (`v*`) is pushed

Untagged builds on `master` produce development versions like `0.10.1.dev20260325+abcd1234`. Tagged releases produce clean version numbers (see below).

## Creating a Release with a Clean Version Number

1. Update the version in `CMakeLists.txt` if needed:
   ```cmake
   set(MADNESS_MAJOR_VERSION 0)
   set(MADNESS_MINOR_VERSION 10)
   set(MADNESS_MICRO_VERSION 2)
   ```

2. Commit the version bump:
   ```
   git commit -am "Bump version to 0.10.2"
   ```

3. Create and push a tag:
   ```
   git tag v0.10.2
   git push origin master --tags
   ```

The tag triggers the conda deploy workflow directly (without waiting for CI), producing a package with version `0.10.2`.

## Setup Requirements

The workflow requires an **Anaconda.org API token** stored as a GitHub repository secret:

1. Go to [anaconda.org](https://anaconda.org) and log in to the `m-a-d-n-e-s-s` organization
2. Navigate to Settings > Access and generate an API token with upload permissions
3. In the GitHub repository, go to Settings > Secrets and variables > Actions
4. Add a new secret named `ANACONDA_TOKEN` with the token value

## Manual Build (Legacy)

The `Makefile.in` and `recipe/` directory contain the original manual workflow. To use it:

1. Configure MADNESS with CMake (this generates `admin/conda/Makefile` from `Makefile.in`)
2. In the build directory's `admin/conda/`, run `make`
3. Upload manually:
   ```
   anaconda login
   anaconda upload /path/to/madness-*.tar.bz2
   ```
