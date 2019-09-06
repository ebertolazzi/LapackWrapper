submodule compilation
=====================

**LapackWrapper**

Prerequisite
------------

~~~~
cd third_parties
~~~~

**Linux**

~~~~
rake install_linux
~~~~

**OSX**

~~~
rake install_osx
~~~

**Windows**

* Visual Studio 2017 64 bit

~~~~
rake install_win[2017,x64]
~~~~

* Visual Studio 2017 32 bit

~~~~
rake install_win[2017,x86]
~~~~

COMPILE
=======

**On linux**

~~~~
make OPENBLAS=1 config
make clean
make
make install_local
~~~~

or using rake

~~~~
rake build_linux
~~~~

**On windows**

using MINGW on a bash shell

~~~~
make OPENBLAS=1 config
make clean
make
make install_local
~~~~

or using Visual Studio

~~~~
rake build_win[2017,x64]
rake build_win[2017,x86]
~~~~

**On OSX use**

~~~~
make ACCELERATE=1 config
make clean
make
make install_local
~~~~

**the library and header file are in the directory `lib`**

~~~~
lib/lib  (static)
   /dll  (dll windows import library)
         (unix shared library)
   /bin  (windows dynamic lib [the dll])
   /headers
~~~~