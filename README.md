# mesh_convert
Adapting meshes to generate anisotropic elements

## Usage

```shell
./mesh_read netgen.vol netgen.vtu metric_file.mtr
```

The file `netgen.vtu` represents the input mesh,
and can then be opened in Paraview.
The directory `adapting` contains a visualization
of adaptation.
Open `adapting/triangles/steps.pvd` in Paraview
to see the mesh during adaptation.
This series will show a `quality` field on triangles.
To see edge lengths, open `adapting/edges/steps.pvd`
in Paraview and look at the `length` field on
"cells" (edges in this series).

All these files should have a `metric` field on
vertices as well.
