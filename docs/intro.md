# dualGPy
Python library to generate the dual graph of a generic mesh.
Starting from the idea of the [meshio](https://github.com/nschloe/meshio) library we developed a tool to read **unstructured meshes**, and give back the [CSR representation](https://en.wikipedia.org/wiki/Sparse_matrix) of the Adjacency Matrix, relative to the dual Mesh. The idea is to use the ASAP (As Simple As Possible) rule.


![Profile Image](https://github.com/albiremo/dualGPy/blob/main/profile.png)

## Mesh class
The `Mesh` class has been developed as an abstact class to be specified for the 2D and 3D (WIP):
- Volume of the cells (Area in 2D);
- Area of the boundary faces (length of the segments in 2D);
- list of cells (starting from the `meshio` format we retrive the cells points);
- list of points composing the mesh.

## Graph class

Built with the same philosophy of `Mesh` class it is able to operate and to build the CSR representation of the Adjacency Matrix.

## Geometry class

A hierarchy of geometrical object to treat the 2D and the 3D meshes. 

## Utils 

Useful utils for the pre and post processing (mesh management, like passing from a 3D to a 2D mesh and output printing)

## Tests

The tests are being run through the [Pytest](https://docs.pytest.org/en/7.1.x/contents.html) suite. If you want to add a test follow the format proposed.

## Examples

to look at an example have a look to the `test.py` script.`

## Planned improvements

To have an idea of the planned improvements, have a look at the issue tracker.

