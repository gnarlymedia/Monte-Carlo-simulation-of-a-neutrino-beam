# Monte-Carlo simulation of a neutrino beam
Monte Carlo methods are computational methods relying on repeated random sampling to obtain numerical results. They can be used to computationally simulate systems which contain an element of chance - an event may or may not occur, or the time or place at which it occurs may follow a probabilistic distribution.

Here we look at a system which produces a neutrino beam through meson decay. A beam consisting of  -mesons and K-mesons (pions and kaons) is allowed to enter a decay tunnel. Here as the beam propagates, it undergoes meson decay into a muon and a neutrino.

## Build and run
Make sure you have cpgplot installed, the C-callable version of pgplot <http://www.astro.caltech.edu/~tjp/pgplot/>

### Command line
Build:
```
make proj_3
```

Then run:
```
./proj_3
```

### Within CLion, JetBrains' IDE for C projects
The `CMakeLists.txt` file should give you everything you need to specify build targets and required libraries, allowing you to build and run within the IDE.