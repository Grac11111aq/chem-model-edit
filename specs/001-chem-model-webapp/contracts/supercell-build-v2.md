# Contract: Supercell Build v2 (Editor v2)

## Endpoints

- POST /structures
- POST /structures/import
- POST /supercell/build

## POST /structures (register Structure)

Request:
```json
{
  "structure": {
    "atoms": [{ "symbol": "O", "x": 0, "y": 0, "z": 0 }],
    "lattice": {
      "a": { "x": 5, "y": 0, "z": 0 },
      "b": { "x": 0, "y": 5, "z": 0 },
      "c": { "x": 0, "y": 0, "z": 5 }
    }
  },
  "source": "editor-v2"
}
```

Response:
```json
{
  "structureId": "9f1c2d7d0c9a4a6c8c5c34f3b8b5b6e3",
  "source": "editor-v2"
}
```

## POST /structures/import (QE import)

Request:
```json
{
  "content": "<QE .in>"
}
```

Response:
```json
{
  "structureId": "a6b5c4d3e2f1a6b5c4d3e2f1a6b5c4d3",
  "structure": {
    "atoms": [{ "symbol": "O", "x": 0, "y": 0, "z": 0 }],
    "lattice": {
      "a": { "x": 5, "y": 0, "z": 0 },
      "b": { "x": 0, "y": 5, "z": 0 },
      "c": { "x": 0, "y": 0, "z": 5 }
    }
  },
  "source": "ase"
}
```

## POST /supercell/build

Request:
```json
{
  "baseStructureId": "11111111111111111111111111111111",
  "grid": {
    "rows": 2,
    "cols": 3,
    "tiles": [
      ["11111111111111111111111111111111", "22222222222222222222222222222222", "11111111111111111111111111111111"],
      ["33333333333333333333333333333333", "11111111111111111111111111111111", "22222222222222222222222222222222"]
    ],
    "axis": { "row": "b", "col": "a" }
  },
  "options": {
    "checkOverlap": true,
    "overlapTolerance": 0.2,
    "validateLattice": "none"
  },
  "output": {
    "includeStructure": false
  }
}
```

Response:
```json
{
  "structureId": "44444444444444444444444444444444",
  "meta": {
    "rows": 2,
    "cols": 3,
    "tileCount": 6,
    "overlapCount": 1,
    "baseStructureId": "11111111111111111111111111111111",
    "structureIdsUsed": [
      "11111111111111111111111111111111",
      "22222222222222222222222222222222",
      "33333333333333333333333333333333"
    ]
  }
}
```

### Field definitions (supercell/build)

- grid.axis: map grid rows/cols to lattice axes (`row` and `col` must be different).
- options.validateLattice:
  - `none` (default): skip lattice validation.
  - `warn`: proceed even if lattices mismatch; server may log a warning.
  - `error`: reject the request if lattices mismatch.
- options.checkOverlap: when true, the build still succeeds and `overlapCount` is returned.

### Error responses

Invalid grid dimensions (400):
```json
{ "detail": "tiles[1] has 2 cols, expected 3" }
```

Unknown structureId (404):
```json
{ "detail": "Structure not found" }
```

Missing base lattice (400):
```json
{ "detail": "base structure has no lattice" }
```

## Notes

- Unit: angstrom (Cartesian coordinates and lattice vectors).
- Tiles are translated based on the lattice of baseStructureId.
- If grid.axis is omitted, row->b and col->a are used by default.
- tiles must be rectangular and match rows/cols.
- structureId references the backend StructureRegistry (ASE Atoms).
- When output.includeStructure=true, the response includes the Structure payload.
- /supercell/build always registers the output structure and returns structureId.
