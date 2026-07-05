# Conda Packaging for MADNESS

Conda packages are built and uploaded to [anaconda.org/m-a-d-n-e-s-s](https://anaconda.org/m-a-d-n-e-s-s) automatically via GitHub Actions.

## Available Packages

| Package | Contents | MPI | Install command |
|---------|----------|-----|-----------------|
| `madness` | C++ libraries, headers, CMake config, apps, data | no | `conda install -c m-a-d-n-e-s-s madness` |
| `madness-mpi` | same, with MPI | Open MPI (from conda-forge) | `conda install -c m-a-d-n-e-s-s -c conda-forge madness-mpi` |
| `pymadness` | Python bindings (depends on `madness`) | no | `conda install -c m-a-d-n-e-s-s -c conda-forge pymadness` |
| `pymadness-mpi` | Python bindings (depends on `madness-mpi`) | Open MPI (from conda-forge) | `conda install -c m-a-d-n-e-s-s -c conda-forge pymadness-mpi` |

All packages are built for Linux and macOS.

The **no-MPI** variant is the simplest option for single-node workstation use.

The **MPI variant** uses OpenMPI from conda-forge and is suitable for workstations and small clusters. Note that the conda-provided OpenMPI uses TCP networking and will not take advantage of high-performance interconnects (InfiniBand, Cray Aries, etc.).

The **`pymadness` packages** carry the Python bindings (see `src/pymadness/`). Each declares an exact dependency on the matching `madness` build, so `conda install pymadness` pulls in `madness` automatically — there is no need to install them separately, and the two are guaranteed to be the same build. The bindings statically embed the MADNESS they were built against; to use a locally compiled MADNESS from Python, build the bindings from source (`cmake -DENABLE_PYTHON=ON`) instead of installing the conda package.

**For HPC clusters**: expert users who need optimal performance should build MADNESS from source to take full advantage of their hardware (native MPI, vendor-tuned BLAS/LAPACK, InfiniBand, etc.). See the main project README for build instructions.

## Conda Recipe

The conda recipe lives in `/conda-recipe/` and consists of:

- **`meta.yaml`** — A **multi-output** recipe: MADNESS is compiled once and the install tree is split into two packages, `madness[-mpi]` (C++ libraries, headers, CMake config, apps, data) and `pymadness[-mpi]` (the Python package + `_pymadness` extension). Uses Jinja templating with `GIT_DESCRIBE_TAG` for versioning and `MPI_VARIANT` to control MPI support; Jinja conditionals toggle the MPI dependencies and the `-mpi` suffix. The `pymadness` output pins its `madness` counterpart exactly via `pin_subpackage(..., exact=True)`.
- **`build.sh`** — Build script invoked by `conda build`. Runs the full CMake build inside the conda environment using conda-forge compilers (with `-DENABLE_PYTHON=ON` and position-independent code so the static libraries can be linked into the extension), relocates the Python package into `site-packages`, and installs activation scripts that set `MRA_DATA_DIR` and `MRA_CHEMDATA_DIR`.

> **Note:** `/conda-recipe/` is the single source of truth for conda packaging and is what the GitHub Actions workflow builds. The templates under `admin/conda/recipe/` (`*.in`) belong to the older manual workflow described under *Manual Build (Legacy)* below and are not used by the automated build.

## Automated Workflow

The workflow is defined in `/.github/workflows/conda-deploy.yml`. It runs automatically when:
- The CI tests pass on `master` (triggered after the "Linux/MacOS Build" workflow succeeds)
- A version tag (`v*`) is pushed
- Manually via `workflow_dispatch`

For each trigger, four build jobs run (2 platforms x 2 MPI variants), and each job emits two conda packages — `madness[-mpi]` and `pymadness[-mpi]` — for eight uploaded packages in total. Build artifacts are also saved as GitHub workflow artifacts.

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
