---
tags:
  - language
  - cheatsheet
---


Cython can be used both with C/C++ syntax and pure Python syntax. It s a prefered extension to Pyrex.

Code is written in a `.pyx` file and compiled with a `setp.py` file. Building is done with:

```bash
python3 setup.py build_ext --inplace
```

https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html
