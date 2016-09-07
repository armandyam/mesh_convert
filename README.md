# mesh_convert
Adapting meshes to generate anisotropic elements

## Usage

```shell
./omega_h netgen.vol metric_file.mtr adapted.vol
Options:
  --adapt_log <vtk_output_directory>
  --axes <output.vtu>
  --limit <scalar (1.0 good default)>
```

The input mesh and metric are specified first;
the third argument is the file to write the adapted mesh to.

If the optional `--adapt_log` argument is specified as
"adapting", the directory `adapting` contains a visualization
of adaptation.
Open `adapting/triangles/steps.pvd` in Paraview
to see the mesh during adaptation.
This series will show a `quality` field on triangles.
To see edge lengths, open `adapting/edges/steps.pvd`
in Paraview and look at the `length` field on
"cells" (edges in this series).
All these files should have a `metric` field on
vertices as well.

If the optional `--axes` argument is specified as `axes.vtu`,
that file will be created with vector fields `axis_0` and
`axis_1` that can be given to a Paraview Glyph filter to render
major and minor axes of the metric field.

If the optional `--limit` argument is specified with value 1.0,
the rate of change of the metric *over metric distance* will
be limited to this rate.
This is useful if you have a few problematic metric values which
are requesting elements larger than can possibly be satisfied.
Desired lengths will be reduced (never increased) as necessary
to limit the metric field to something more satisfiable.
