# PCGP

PCGP is a multithreaded C/C++ program for counting the bisection of circulant graphs using Kernighan-Lin algorithm.

This is a singlethreaded implementation.

# Usage

## Windows

### Build

Compile all .c files with any C99 compiler. No dependencies.

### Run

Edit parameters N, K, C in [pcgp.bat](/pcgp.bat) script and execute it.

## Linux

Run the following command:

```bash
make N=<n> K=<k> C=<c>
```

Change `<n>`, `<k>`, `<c>` to the desired parameter values.

This command will build the executable if it is not built yet and run it with the defined parameters.

