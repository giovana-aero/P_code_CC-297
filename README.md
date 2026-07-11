# \<\< CC-297 - Elementos de mecânica dos fluidos computacional - Project code \>\>

Folder structure:
- P01: Biconvex airfoil, transonic small disturbance potential solution
- P02: Elliptical type O mesh over a biconvex airfoil
- P03: Full potential flow solution over a biconvex airfoil
- exercise_laplace2d: solution of the two-dimensional Laplace equation. I wrote this as practice before working on P01, intending to reuse the functions in the projects. Unfortunately, the code was only partially reusable.
- post: generic post-processing scripts. Specialized scripts for each project are saved in their respective folders.
- src: source code of the libraries.
- include: header files of the libraries.

Explanations regarding the configuration parameters are available in the `include` files. To compile a given case, compile the main script along with the source, such as:

`$ gcc biconvex_airfoil.c ../src/* -lm    # from inside the P01 folder`

Where the `-lm` flag is required by the `math.h` library.

There are some MATLAB functions I have not included here yet.


---

Giovana W. Fernandes, MSc\
Instituto Tecnológico de Aeronáutica\
São josé dos Campos, 2026
