/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Sparse_tool  : Simple Class for timing code                             |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : Sparse_tool_Extra.hh                                     |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#pragma once
#ifndef SPARSETOOL_EXTRA_HH
#define SPARSETOOL_EXTRA_HH

#include "sparse_tool.hh"

namespace Sparse_tool {

  /*
  //   ____
  //  / ___| _ __  _   _
  //  \___ \| '_ \| | | |
  //   ___) | |_) | |_| |
  //  |____/| .__/ \__, |
  //        |_|    |___/
  */

  //!
  //! Plot a sihouette of the sparse pattern to a file
  //! \param fname    file name
  //! \param sp       sparse pattern to plot
  //! \param xsize    dimension of the plot in `cm`
  //! \param rowLines if not `nullptr` is a pointer to a `vector` object
  //!                 with the list of the horizontal lines to draw
  //! \param colLines if not `nullptr` is a pointer to a `vector` object
  //!                 with the list of the vertical lines to draw
  //!
  template <typename MAT>
  void
  Spy(
    string_view             fname,
    SparseBase<MAT> const & sp,
    double          const   xsize,
    vector<integer> const * rowLines = nullptr,
    vector<integer> const * colLines = nullptr
  ) {
     std::ofstream file( fname.data() );
     Spy( file, sp, xsize, rowLines, colLines );
     file.close();
  }

  //!
  //! Plot a sihouette of the sparse pattern to a file
  //! \param stream   stream object file where the figure is saved
  //! \param sp       sparse pattern to plot
  //! \param xsize    dimension of the plot in `cm`
  //! \param rowLines if not `nullptr` is a pointer to a `vector` object
  //!                 with the list of the horizontal lines to draw
  //! \param colLines if not `nullptr` is a pointer to a `vector` object
  //!                 with the list of the vertical lines to draw
  //!
  template <typename MAT>
  void
  Spy(
    ostream_type          & stream,
    SparseBase<MAT> const & sp,
    double          const   xsize,
    vector<integer> const * rowLines = nullptr,
    vector<integer> const * colLines = nullptr
  ) {

    /*
    // based on:
    //
    // PSPLTM - PostScript PLoTer of a (sparse) Matrix
    // By Loris Renggli (renggli@masg1.epfl.ch) and Youcef Saad
    */

    // ----------------------------------------------------------------------
    integer nc = sp.ncols() + 1;
    integer nr = sp.nrows() + 1;

    double ysize = (xsize*nr) / nc;

    // units (cm or in) to dot conversion factor and paper size
    double u2dot  = 72.0/2.54;
    double paperx = 21.0;
    double papery = 29.7;

    // left and right margins (drawing is centered)
    double leftRightMargin = (paperx-xsize)/2.0;
    double bottomTopMargin = (papery-ysize)/2.0;

    // scaling factor
    //double scalingFactor = std::min(xsize/nc,ysize/nr);

    // matrix frame line witdh
    double sepLaneWidth   = 0.1;

    // matrix frame line witdh
    double frameLaneWidth = std::max(sepLaneWidth,0.001*nc);

    // almost exact bounding box
    double xl = leftRightMargin;
    double xr = xl+xsize;
    double yb = bottomTopMargin;
    double yt = yb + ysize;

    // add some room to bounding box and convert to dots
    double delt = 5.0;
    xl = xl * u2dot - delt;
    xr = xr * u2dot + delt;
    yb = yb * u2dot - delt;
    yt = yt * u2dot + delt;

    // begin of output
    fmt::print( stream,
      "%!PS-Adobe-3.0 EPSF-3.0\n"
      "%%Creator: Sparse_tool by Enrico Bertolazzi\n"
      "%%BoundingBox: {} {} {} {}\n"
      "%%EndComments\n"
      "/cm {{72 mul 2.54 div}} def\n"
      "gsave\n"
      "{} cm {} cm translate\n"
      "{} cm {} div dup scale\n"
      // draw a frame around the matrix
      "{} setlinewidth\n"
      "newpath\n"
      "{} {} moveto\n"
      "{} {} lineto\n"
      "{} {} lineto\n"
      "{} {} lineto\n"
      "closepath stroke\n",
      round(xl), round(yb), round(xr), round(yt),
      leftRightMargin, bottomTopMargin,
      xsize, nc,
      frameLaneWidth,
      0.5-frameLaneWidth,    0.5-frameLaneWidth,
      nc-0.5+frameLaneWidth, 0.5-frameLaneWidth,
      nc-0.5+frameLaneWidth, nr-0.5+frameLaneWidth,
      0.5-frameLaneWidth,    nr-0.5+frameLaneWidth
    );

    // drawing the separation lines (if required)
    vector<integer>::const_iterator it;
    if ( rowLines != nullptr ) {
      fmt::print( stream, "{} setlinewidth\n", sepLaneWidth );
      for ( it = rowLines -> begin(); it != rowLines -> end(); ++it ) {
        double y = sp.nrows() - (*it) + 0.5;
        fmt::print( stream, "0.5 {} moveto {} {} lineto stroke\n", y, nc-0.5, y );
      }
    }
    if ( colLines != nullptr ) {
      fmt::print( stream, "{} setlinewidth\n", sepLaneWidth );
      for ( it = colLines -> begin(); it != colLines -> end(); ++it ) {
        double x = (*it)+0.5;
        fmt::print( stream, "{} 0.5 moveto {} {} lineto stroke\n", x, x, nr-0.5 );
      }
    }

    // ----------- plotting loop
    fmt::print( stream,
      "1 1 translate\n"
      "0.8 setlinewidth\n"
      "/p {{ moveto 0 -.40 rmoveto 0 .80 rlineto stroke }} def\n"
    );

    for ( sp.Begin(); sp.End(); sp.Next() )
      fmt::print( stream, "{} {} p\n", sp.column(), (sp.nrows()-sp.row()-1) );

    // -----------------------------------------------------------------------
    stream << "showpage\n";
  }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Sparse_tool_load {
  using ::Sparse_tool::Spy;
}
#endif

#endif
