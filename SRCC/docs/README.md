# API Documentation of C++ source code

The API documentation of the C++ source code is generated from the source code using
`doxygen`.

## Build

In order to build the API documentation, you need to
[install doxygen](https://www.doxygen.nl/download.html).
Once installed from the root of the repository:

```bash
cd SRCC
doxygen docs/Doxyfile
```

This will build the documentation in the `SRCC/docs/html` directory.

## View

Simply open the `SRCC/docs/html/index.html` file in your browser or consult the
latest [online documentation](https://payoto.github.io/rsvs3d/index.html).

## Dependencies

The doxygen theme used by this repository is
[Doxygen Awesome](https://github.com/jothepro/doxygen-awesome-css).
