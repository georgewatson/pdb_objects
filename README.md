# PDB Objects

Set of classes for object-oriented processing of records from PDB files

© 2019  George D. Watson, University of York
([https://georgewatson.me](https://georgewatson.me))

Available under an MIT license. See the LICENSE file.

This package is described in
[a post on my blog](https://georgewatson.me/blog/science/2019/06/04/object-oriented-processing-of-pdb-files-in-python/).

Supports the following record types:
* `ATOM`
* `HETATM`
* `TER`
* `HELIX`
* `SHEET`

No other record types are currently implemented, but can be implemented upon
request.

Exposes the following classes:
* `PDBRecord` (should not normally be used directly, except to implement
  another record type)
* `Residue`
* `Atom` (for ATOM and HETATM records)
* `Helix`
* `Sheet`
* `Ter`

Exposes the following public functions:
* `read_atom`
* `read_helix`
* `read_sheet`
* `read_ter`
* `read_record`
* `read_pdb`

See class, function, and method docstrings for more information.

## Installation

Available on PyPI.
Use `pip3 install pdb-objects` to install,
then put `import pdb_objects` at the top of your script.

Alternatively, clone this repository.
