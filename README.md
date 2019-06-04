# PDB Objects

Set of classes for object-oriented processing of records from PDB files

© 2019  George D. Watson, University of York
[https://georgewatson.me](https://georgewatson.me)

Available under an MIT license. See the LICENSE file.

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
