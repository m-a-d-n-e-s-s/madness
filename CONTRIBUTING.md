# Contributing to MADNESS

## License of contributions

MADNESS is distributed under the MIT License (see `LICENSE`).
Contributions are accepted under the same license ("inbound = outbound"):
by submitting a pull request you agree that your contribution is licensed
under the MIT License.

You keep the copyright to your contribution. If your employer owns the
copyright to work you do, make sure you are allowed to contribute it under
these terms before you submit it.

## Sign-off

Every commit must carry a `Signed-off-by:` line certifying the
[Developer Certificate of Origin](https://developercertificate.org/):

```
git commit -s
```

The sign-off states that you wrote the change, or otherwise have the right
to submit it under the project's license. It is the project's record of
where each contribution came from.

## Copyright headers

Each source file starts with a short header:

```cpp
/*
  This file is part of MADNESS.

  SPDX-FileCopyrightText: 2007, 2010 Oak Ridge National Laboratory
  SPDX-FileCopyrightText: The MADNESS developers
  SPDX-License-Identifier: MIT
*/
```

- **New files:** use the header above with only the
  `The MADNESS developers` copyright line, followed by
  `SPDX-License-Identifier: MIT`.
- **Editing a file:** you do not need to change the header or update years.
  Git history records who changed what.
- **Institutional notices:** if you write a file, or a substantial part of
  one, and your institution requires its own notice, add an
  `SPDX-FileCopyrightText: <year> <institution>` line above the
  `The MADNESS developers` line.
- **Never remove** an existing copyright line. Only the holder it names can
  agree to its removal.
- **Third-party code** (for example under `src/madness/external`) keeps its
  original header and license.

## Contributors list

`CONTRIBUTORS.md` acknowledges the people who contributed to MADNESS while it
was distributed under the GPL. Contributors since the relicensing are
recorded in the git history.
