SparseTool: a Sparse Matrix Manager
===================================

Preface
-------

``SparseTool`` is a collection of simple, fast, and efficient classes
for manipulating large vectors and large sparse matrices.

The Basic Vector and Matrix Classes
-----------------------------------

The library consist of the following templated classes: - ``Vector<T>``
which define a **dense** large column vector. - ``SparsePattern`` which
define a **pattern** of the nonzero elements which can be used to
construct a sparse matrix. - ``TridMatrix<T>`` which implements a
**tridiagonal** matrix. - ``CCoorMatrix<T>`` which implements a sparse
``Compressed Coordinate Storage`` matrix. - ``CRowMatrix<T>`` which
implements a sparse ``Compressed Rows Storage`` matrix. -
``CColMatrix<T>`` which implements a sparse
``Compressed Columns Storage`` matrix.

Those classes allow vectors and matrices to be formally treated in
software implementations as mathematical objects in arithmetic
expressions. For example if ``A`` is a sparse matrix and ``b`` is a
``Vector<T>`` than ``A*b`` means the matrix-vector product.

Loading the library
~~~~~~~~~~~~~~~~~~~

To use the library you must include it by the following piece of code:

.. code:: cpp

   #include "SparseTool.hh"
   using namespace SparseToolLoad;

line **2** is recommended to avoid ``SparseTool::`` prefix for the
library call.

Some Macros for customization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To acivate some control when ``SparseTool``\ is used in developing

.. code:: cpp

   #define SPARSETOOL_DEBUG

Variable size full vector class
-------------------------------

This class extend the STL vector class by adding some math operation and
interaction with sparse matrix classes.

Usage
~~~~~

A ``Vector<T>`` is defined by specifying the type ``T`` and optionally
its size. For example,

.. code:: cpp

     Vector<double> b, c(100);

defines ``b`` as a vector of **double** of size **0** while ``c`` is a
vector of **double** of size **100**. You can change the size of the
vector by the methods ``resize`` as follows:

.. code:: cpp

   Vector<double> d;
   d.resize(200);

so that ``d`` is a ``Vector<double>`` of size **200**. There are many
methods associate to a ``Vector<T>`` in the following paragraph they are
listed.

Constructors
~~~~~~~~~~~~

.. code:: cpp

   Vector<T> v;
   Vector<T> v(dim);
   Vector<T> w(v);

On line **1** construct the ``Vector<T> v`` of size **0**. On line **2**
construct the ``Vector<T> v`` of of size ``dim``. On line **3**
construct the ``Vector<T> w`` as a copy of ``Vector<T> v``.

Indexing
~~~~~~~~

Vector instances are indexed as one-dimensional C arrays, and the index
numbering follows the standard C convention, starting from zero.

Let us define ``v`` as ``Vector<T>``

.. code:: cpp

   Vector<T> v;

Then, ``v[i]`` returns a reference to the ``T&-type`` ``i-th`` element
of ``v``;

Changing Dimension
~~~~~~~~~~~~~~~~~~

It is possible to change the size of a ``Vector<T>`` object. For example

.. code:: cpp

   Vector<double> d;
   d.resize(200);

the ``Vector<T> d`` has size **200**. The method ``size()`` return the
actual size of the ``Vector<T>``. For example defining

.. code:: cpp

   Vector<double> d(123);

the method ``d.size()`` return **123**.

Initialization
~~~~~~~~~~~~~~

It is possible to initialize all the components of a ``Vector<T>`` to a
value, for example

.. code:: cpp

   Vector<float> v(100);
   v = 3.14;

is equivalent to

.. code:: cpp

   Vector<float> v(100);
   for ( int i = 0; i < v.size(); ++i ) v[i] = 3.14;

although is done more efficiently by the library.

Assignment
~~~~~~~~~~

It is possible to copy the contents of a ``Vector<T>`` to another one as
the following example show

.. code:: cpp

   Vector<double> a, b, c;

   a.resize(100);
   b.resize(200);
   c.resize(150);

   c = 3;
   b = c;
   a = c;

the meaning of line **7** should be clear. Lines **8** and **9** are
equivalent to

.. code:: cpp

   int i;
   for ( i = 0; i < min(b.size(), c.size()); ++i ) b[i] = c[i];
   for ( i = 0; i < min(a.size(), c.size()); ++i ) a[i] = c[i];

where you can notice that only the values that can be stored are
assigned.

It is possible to initialize many vectors same value as in the following
expressions

.. code:: cpp

   Vector<double> a, b, c;
   a.resize(100);
   b.resize(200);
   c.resize(150);
   a = b = c = 3;

but take attention because line **5** is not equivalent to

.. code:: cpp

   a = 3;
   b = 3;
   c = 3;

in fact line **5** is equivalent to

.. code:: cpp

   c = 3;
   b = c;
   a = b;

so that we have

-  the ``Vector<T> c`` is initialized with **all** its **150** elements
   set to **5**.
-  the ``Vector<T> b`` is initialized with **only** its first **150**
   elements set to **1** while the remaining are undefined
-  the ``Vector<T> a`` is initialized with **all** its first **100**
   elements set to **1**.

Arithmetic Operators on ``Vector<T>``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A set of usual arithmetic operators are explicitly defined on
vector-type data. If not otherwise specified, the operators extend the
corresponding scalar operation in a **component-wise** fashion. Hence,
for vectors with size ``dim``, the component index ``i`` in all the
following expressions is supposed to run through **0** to ``dim-1``.

Let us define the three double precision vectors ``a``, ``b``, and
``c``, that we shall use in all the following examples

.. code:: cpp

   int const dim = 100;
   Vector<double> a(dim), b(dim), c(dim);

The arithmetic operators defined on vectors are given in the following
sections.

**Scalar-Vector** internal operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

========== =============================================
Command    Equivalence
========== =============================================
``a += 2`` ``for ( i=0; i < a.size(); ++i ) a[i] += 2;``
``a -= 2`` ``for ( i=0; i < a.size(); ++i ) a[i] -= 2;``
``a *= 2`` ``for ( i=0; i < a.size(); ++i ) a[i] *= 2;``
``a /= 2`` ``for ( i=0; i < a.size(); ++i ) a[i] /= 2;``
========== =============================================

**Scalar-Vector** operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

============= =============================================
Command       Equivalence
============= =============================================
\             ``sz = min( a.size(), b.size() );``
``a = b + 2`` ``for ( i=0; i < sz; ++i ) a[i] = b[i] + 2;``
``a = 3 + b`` ``for ( i=0; i < sz; ++i ) a[i] = 3 + b[i];``
``a = b - 2`` ``for ( i=0; i < sz; ++i ) a[i] = b[i] - 2;``
``a = 3 - b`` ``for ( i=0; i < sz; ++i ) a[i] = 3 - b[i];``
``a = b * 2`` ``for ( i=0; i < sz; ++i ) a[i] = b[i] * 2;``
``a = 3 * b`` ``for ( i=0; i < sz; ++i ) a[i] = 3 * b[i];``
``a = b / 2`` ``for ( i=0; i < sz; ++i ) a[i] = b[i] / 2;``
``a = 3 / b`` ``for ( i=0; i < sz; ++i ) a[i] = 3 / b[i];``
============= =============================================

**Vector-Vector** internal operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

========== ==========================================
Command    Equivalence
========== ==========================================
\          ``sz = min( a.size(), b.size() );``
``a += b`` ``for ( i=0; i < sz; ++i ) a[i] += b[i];``
``a -= b`` ``for ( i=0; i < sz; ++i ) a[i] -= b[i];``
``a *= b`` ``for ( i=0; i < sz; ++i ) a[i] *= b[i];``
``a /= b`` ``for ( i=0; i < sz; ++i ) a[i] /= b[i];``
========== ==========================================

**Vector-Vector**
~~~~~~~~~~~~~~~~~

============= ================================================
Command       Equivalence
============= ================================================
\             ``sz = min( a.size(), b.size() );``
``b = +a``    ``for ( i=0; i < sz; ++i ) b[i] = +a[i];``
``b = -a``    ``for ( i=0; i < sz; ++i ) b[i] = -a[i];``
\             ``sz = min( sz, c.size() );``
``c = a + b`` ``for ( i=0; i < sz; ++i ) c[i] = a[i] + b[i];``
``c = a - b`` ``for ( i=0; i < sz; ++i ) c[i] = a[i] - b[i];``
``c = a * b`` ``for ( i=0; i < sz; ++i ) c[i] = a[i] * b[i];``
``c = a / b`` ``for ( i=0; i < sz; ++i ) c[i] = a[i] / b[i];``
============= ================================================

Function of ``Vector<T>``
~~~~~~~~~~~~~~~~~~~~~~~~~

Let be ``n = min( a.size(), b.size() )``,

-  ``dot(a,b)``   :math:`=\displaystyle\sum_{i=0}^{n-1}\overline{a_i}b_i`
-  ``rdot(a,b)``  :math:`=\displaystyle\sum_{i=0}^{n-1} a_i b_i`
-  ``dist(a,b)``  :math:`=\sqrt{\sum_{i=0}^{n-1} \left|a_i - b_i\right|^2}`
-  ``normi(a)``   :math:`= ||a||_\infty =\max\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\}`
-  ``norm1(a)``   :math:`= ||a||_1 =\sum_{i=0}^{n-1}|a_i|`
-  ``norm2(a)``   :math:`= ||a||_{2}=\sqrt{\sum_{i=0}^{n-1} a_i^2}`
-  ``normp(a,p)`` :math:`=||a||_p = \left(\sum_{i=0}^{n-1}|a_i|^p\right)^{1/p}`
-  ``sum(a)``     :math:`=\displaystyle\sum_{i=0}^{n-1} a_i`
-  ``prod(a)``    :math:`=\displaystyle\prod_{i=0}^{n-1} a_i`
-  ``max(a)``     :math:`=\max\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\}`
-  ``min(a)``     :math:`=\min\left\{|a_0|,|a_1|,\ldots,|a_{n-1}|\right\}`

Remapping a piece of memory to a vector
---------------------------------------

This class extend the ``VectorBase`` class by adding the capacity of
remapping a piece o ``Vector`` or a C-pointer

Usage
~~~~~

A ``VectorSlice<T>`` is defined by specifying the type ``T``. For
example,

.. code:: cpp

   VectorSlice<double> b, c;

There are many methods associate to a ``Vector<T>`` in the following
paragraph they are listed.

Slicing
~~~~~~~

.. code:: cpp

   Vector<double> v(100); double w[100];
   b.slice( v, 10, 20 );
   c.slice( w + 5, w + 45 );

On line **1** construct the ``Vector<double>`` ``v`` of size **100** and
the C array ``w`` of **100** elements.

On line **2** remap the components of the vecotor ``v`` from **10** to
**19** to the “vector” ``b``.

On line **3** remap the components of the Array ``w`` from **5** to
**44** to the “vector” ``c``. ``Vector<T> v``.

.. _indexing-1:

Indexing
~~~~~~~~

Vector instances are indexed as one-dimensional C arrays, and the index
numbering follows the standard C convention, starting from zero.

Then, ``b[i]`` returns a reference to the ``T&-type`` ``i-th`` element
of ``v`` and due to slicing the ``b[i]==v[i+10]``. Analogously
``c[i]==w[i+5]``.

The Preconditioner Classes and Iterative algorithms
---------------------------------------------------

The library consists of the following templated classes:

-  ``IdPreconditioner<T>`` which implements the identity preconditioner.
-  ``Dpreconditioner<T>`` which implements the diagonal preconditioner.
-  ``ILDUpreconditioner<T>`` which implement an incomplete LDU preconditioner.

A set of template iterative solvers are available:

-  ``cg`` implementing the cojugate gradient solver
-  ``bicgstab`` implementing the Bi-conjugate stabilized solver of Van
   Der Vorst.
-  ``gmres`` implementing generalized minimal residual of Saad-Shultz

Iterative solvers
~~~~~~~~~~~~~~~~~

Here is an example of the use of the iterative solver:

.. code:: cpp


     CRowMatrix<double> A;
     Vector<double>     b, x;
     ILDUpreconditioner P;
     double             tolerance;
     int                maxIter, iter;
     .
     .
     .
     .
     double residual = cg(A, b, x, P, tolerance, maxIter, iter);

     double residual = bicgstab(A, b, x, P, tolerance, maxIter, iter);

     double residual = gmres(A, b, x, P, tolerance, maxSubIter, maxIter, iter);

In the example

-  ``A`` : is the coefficients matrix;
-  ``b`` : is the known vector;
-  ``x`` : is the vector which will contains the solution;
-  ``P`` : is the preconditioner object class;
-  ``tolerance`` : is the admitted tolerance;
-  ``maxSubIter`` : for ``gmres`` is the maximum number of iteration
   before restarting;
-  ``maxIter`` : is the maximum number of allowable iterations;
-  ``iter`` : is the number of iterations done;
-  ``residual`` : the residual of the approximated solution;

MatrixMarket
------------

A standard way to store and exchange large sparse matrix is to save
matrix accordingly to some **file exchange format**.

For sparse matrices there are two widely used format the Harwell-Boeing
(HB) Sparse Matrix Format and Matrix Market (MM) Sparse Matrix Format.

The HB format is strongly depended on ``FORTRAN`` language and is
difficult to manage in other languages unless using a sophisticated
parser.

The alternative is to use dedicated ``FORTRAN`` routine to manage such a
format.

The design of ``SparseTool``\ is to use a unique header
``SparseTool.hh`` which contains the whole library so that HB format is
not supported.

The MM format is easier to manage so that a simple support is included
in the toolkit by the class ``MatrixMarket``.

In any case there are free software for convert from one format to
another (see e.g. http://bebop.cs.berkeley.edu/).

The following code shows how ``MatrixMarket`` should be used:

.. code:: cpp

   MatrixMarket mm; // define the object mm to manage Matrix Market file
   CCoorMatrix<double> A;   // an empty sparse CCOOR matrix
   SparsePattern sp;        // an empty sparse pattern
   mm.read("hor__131.mtx"); // read matrix and store in mm object
   mm.load_pattern( sp );   // extract the pattern
   mm.load_matrix( A );     // copy loaded matrix in sparse matrix A

The class ``MatrixMarket`` is not a ``template`` class because the value
type of the nonzeros is specified by the format.

The available type in MM format are ``int``, ``double`` and
``complex<double>``. The representation of sparse matrices in MM format
is of type CCOOR so that ``MatrixMarket`` store in this form the loaded
matrix. Full matrices are stored in column major order. Sometimes it can
be useful to access directly the structure of the loaded matrix with the
following methods:

In the library it is not provided a ``write`` method to MM format. This
choice is due to the fact that write a MM file is very easy using the
iterators and the only complication is to write the first few rows of
the header file.
