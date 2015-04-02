R package, aSPU
=================

Wei Pan, Il-Youp Kwak <ikwak@umn.edu>

R/aSPU is an R package for Genetic association testing methods such as aSPU, aSPUw, aSPUpath, GATES, GATE-Simes, HYST etc.



### installation
with `devtools `:
```S
library(devtools)
install_github("ikwak2/aSPU")
```

### Window users

Installation in window machine is a bit difficult. Read instruction. [here](http://cran.r-project.org/doc/manuals/R-admin.html#The-Windows-toolset)


Here are list of things you need to install.

[Rtools](http://cran.us.r-project.org/bin/windows/Rtools/) Downlaod the version that matched with your R.
The installer Rtools*.exe optionally sets the path to components that it installs. So, select the option that add path.

[MikTeX](http://miktex.org/download)

[The Inno Setup installer](http://jrsoftware.org/)

[The command line tools](http://www.cygwin.com/)

[The MinGW-w64 toolchain](http://sourceforge.net/projects/mingw-w64/)


Next difficulty is setting path correctly.

In R [documentation](http://cran.r-project.org/doc/manuals/R-admin.html#The-command-line-tools), an example for a 32-bit build is,
PATH=c:\Rtools\bin;c:\Rtools\gcc-4.6.3\bin;c:\MiKTeX\miktex\bin;
     c:\R\R-3.1\bin\i386;c:\windows;c:\windows\system32

[How to set up path in windows](http://www.computerhope.com/issues/ch000549.htm)

Check the directory you installed the program and set the path.

Hope this work!
