/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : Simple Class for timing code                             |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : SparseTool_Extra.hh                                      |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_EXTRA_HH
#define SPARSETOOL_EXTRA_HH

#include "sparse_tool.hh"

namespace SparseTool {

  /*
  //   ____
  //  / ___| _ __  _   _
  //  \___ \| '_ \| | | |
  //   ___) | |_) | |_| |
  //  |____/| .__/ \__, |
  //        |_|    |___/
  */

  /*!
   *  Plot a sihouette of the sparse pattern to a file
   *  \param fname    file name
   *  \param sp       sparse pattern to plot
   *  \param xsize    dimension of the plot in \c cm
   *  \param rowLines if not \c nullptr is a pointer to a \c vector object
   *                  with the list of the horizontal lines to draw
   *  \param colLines if not \c nullptr is a pointer to a \c vector object
   *                  with the list of the vertical lines to draw
   */

  template <typename MAT>
  void
  Spy(
    std::string       const & fname,
    SparseBase<MAT>   const & sp,
    double            const   xsize,
    vector<indexType> const * rowLines = nullptr,
    vector<indexType> const * colLines = nullptr
  ) {
     std::ofstream file( fname.c_str() );
     Spy( file, sp, xsize, rowLines, colLines );
     file . close();
  }

  /*!
   *  Plot a sihouette of the sparse pattern to a file
   *  \param stream   stream object file where the figure is saved
   *  \param sp       sparse pattern to plot
   *  \param xsize    dimension of the plot in \c cm
   *  \param rowLines if not \c nullptr is a pointer to a \c vector object
   *                  with the list of the horizontal lines to draw
   *  \param colLines if not \c nullptr is a pointer to a \c vector object
   *                  with the list of the vertical lines to draw
   */
  template <typename MAT>
  void
  Spy(
    ostream                 & stream,
    SparseBase<MAT>   const & sp,
    double            const   xsize,
    vector<indexType> const * rowLines = nullptr,
    vector<indexType> const * colLines = nullptr
  ) {

    /*
    // based on:
    //
    // PSPLTM - PostScript PLoTer of a (sparse) Matrix
    // By Loris Renggli (renggli@masg1.epfl.ch) and Youcef Saad 
    */

    // ----------------------------------------------------------------------
    indexType nc = sp.numCols() + 1;
    indexType nr = sp.numRows() + 1;

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
    stream
      << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%Creator: SparseTool by Enrico Bertolazzi\n"
      << "%%BoundingBox: " << round(xl) << ' '
                           << round(yb) << ' '
                           << round(xr) << ' '
                           << round(yt) << '\n'
      << "%%EndComments\n"
      << "/cm {72 mul 2.54 div} def\n"
      << "gsave\n"
      << leftRightMargin << " cm " << bottomTopMargin << " cm translate\n"
      << xsize << " cm " << nc << " div dup scale\n"
      // draw a frame around the matrix
      << frameLaneWidth << " setlinewidth\n"
      << "newpath\n"
      << 0.5-frameLaneWidth    << ' ' << 0.5-frameLaneWidth    << " moveto\n"
      << nc-0.5+frameLaneWidth << ' ' << 0.5-frameLaneWidth    << " lineto\n"
      << nc-0.5+frameLaneWidth << ' ' << nr-0.5+frameLaneWidth << " lineto\n"
      << 0.5-frameLaneWidth    << ' ' << nr-0.5+frameLaneWidth << " lineto\n"
      << "closepath stroke\n";

    // drawing the separation lines (if required)
    vector<indexType>::const_iterator it;
    if ( rowLines != nullptr ) {
      stream << sepLaneWidth << " setlinewidth\n";
      for ( it = rowLines -> begin(); it != rowLines -> end(); ++it ) {
        double y = sp.numRows() - (*it) + 0.5;
        stream << "0.5 " << y << " moveto " << nc-0.5 << ' ' << y << " lineto stroke\n";
      }
    }
    if ( colLines != nullptr ) {
      stream << sepLaneWidth << " setlinewidth\n";
      for ( it = colLines -> begin(); it != colLines -> end(); ++it ) {
        double x = (*it)+0.5;
        stream << x << " 0.5 moveto " << x << ' ' << nr-0.5 << " lineto stroke\n";
      }
    }

    // ----------- plotting loop
    stream
      << "1 1 translate\n"
      << "0.8 setlinewidth\n"
      << "/p { moveto 0 -.40 rmoveto 0 .80 rlineto stroke } def\n";

    for ( sp.Begin(); sp.End(); sp.Next() )
      stream << sp.column() << ' ' << (sp.numRows()-sp.row()-1) << " p\n";

    // -----------------------------------------------------------------------
    stream << "showpage\n";
  }

}

namespace SparseToolLoad {
  // utilities
  using ::SparseTool::Spy;
}

#endif
