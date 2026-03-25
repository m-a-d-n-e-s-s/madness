# Conda Packaging for MADNESS

Conda packages are built and uploaded to [anaconda.org/m-a-d-n-e-s-s](https://anaconda.org/m-a-d-n-e-s-s) automatically via GitHub Actions.

## Available Packages

| Package | MPI | Install command |
|---------|-----|-----------------|
| `madness` | no | `conda install -c m-a-d-n-e-s-s madness` |
| `madness-mpi` | OpenMPI (from conda-forge) | `conda install -c m-a-d-n-e-s-s -c conda-forge madness-mpi` |

Both packages are built for Linux and macOS.

The **no-MPI** variant is the simplest option for single-node workstation use.

The **MPI variant** uses OpenMPI from conda-forge and is suitable for workstations and small clusters. Note that the conda-provided OpenMPI uses TCP networking and will not take advantage of high-performance interconnects (InfiniBand, Cray Aries, etc.).

**For HPC clusters**: expert users who need optimal performance should build MADNESS from source to take full advantage of their hardware (native MPI, vendor-tuned BLAS/LAPACK, InfiniBand, etc.). See the main project README for build instructions.

## Automated Workflow

The workflow is defined in `/.github/workflows/conda-deploy.yml`. It runs automatically when:
- The CI tests pass on `master` (triggered after the "Linux/MacOS Build" workflow succeeds)
- A version tag (`v*`) is pushed

For each trigger, four packages are built and uploaded (2 platforms x 2 MPI variants).

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

The tag triggers the conda deploy workflow directly, producing packages with version `0.10.2`.

## Setup Requirements

The workflow requires an **Anaconda.org API token** stored as a GitHub repository secret:

1. Go to [anaconda.org](https://anaconda.org) and log in to the `m-a-d-n-e-s-s` organization
2. Navigate to Settings > Access and generate an API token with upload (write) permissions
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
