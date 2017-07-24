## C++ implementation of pairwise distances in ROG

- implementation uses generator of `combinations_with_repetition` by Siegfried Koepf
http://www.aconnect.de/friends/editions/computer/combinatoricode_e.html

- compilation (`boost` required)
```bash
$ c++ -o rog_lengths rog_lengths.cpp comb_rep_lex.c -ffast-math -march=native -std=gnu++0x -O3
```

- run for specific n and nvar
```bash
$ ./rog_lengths n nvar
$ ./rog_lengths 2 3
```

- read results stored in binary files using `numpy` package

```python
import numpy as np

n = 3
nvar = 2

lens = np.fromfile('lens_ae_{:06d}_{:02d}.bin'.format(n, nvar), dtype=float)
freq = np.fromfile('freq_ae_{:06d}_{:02d}.bin'.format(n, nvar), dtype=float)
```
