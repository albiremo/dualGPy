## HOW TO BUILD THE DOC

if you need to modify anything in the apidoc run

```
sphinx-apidoc -o . ../dualGPy/
```

otherwise to compile how it is run:

```
make html
```

the artifacts will be in the `__build` directory under the name `index.html`.
