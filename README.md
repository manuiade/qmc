# Quine–McCluskey algorithm

Implementation in C of the Quine-McCluskey Algorithm, which performs minimization of Boolean Functions in the form of Sum of Produts or Product of Sums.

<b> That was my first complex algorithm written in C a long time ago, so sorry if the code is really bad to read!! </b>

## Usage

Compile the program:

```
gcc main.c qmc.c -o qmc
```

Launch it:

```
./qmc
```

Specify the number of inputs, the type of minimization you want (Sum of Products or Product of Sums) and the boolean outputs for each possible configuration.