//-----------------------------------------------------------------------------
// File:        PDFweight.cc
// Created:     05-June-2014 - reworking of some old code
// Author:      Harrison B. Prosper      
//-----------------------------------------------------------------------------
#include <iostream>
#include "LHAPDF/LHAPDF.h"
#include "PDFweight.h"
//-----------------------------------------------------------------------------
using namespace std;
//-----------------------------------------------------------------------------
PDFweight::PDFweight(string pdfName, string pdfNameOrig)
  : pdfset_(LHAPDF::PDFSet(pdfName)),
    pdf_(pdfset_.mkPDFs()),
    pdf0_(pdfNameOrig!="" ? LHAPDF::mkPDF(pdfNameOrig, 0) : pdf_[0])
{
  cout << "\t==> using PDF set:      " << pdfset_.name() << endl
       << "\t==> number of PDF sets: " << pdfset_.size() << endl
       << "\t==> original PDF set:   " << pdfNameOrig << endl;
}


PDFweight::~PDFweight() 
{
  if ( pdf_[0] != pdf0_ ) delete pdf0_;
  for(size_t i=0; i < pdf_.size(); i++) 
    delete pdf_[i];
}

double PDFweight::operator()(int id1, int id2,
			     double x1, double x2, 
			     double Q, 
			     int number)
{
  if ( number < 0 ) return -1;
  if ( number > (int)pdfset_.size() ) return -2;

  double pdf1_0, pdf2_0;
  double pdf1_1, pdf2_1;

  // Get central value of PDFs
  pdf1_0 = pdf0_->xfxQ(id1, x1, Q);
  pdf2_0 = pdf0_->xfxQ(id2, x2, Q);

  // Get non-central value of PDFs
  pdf1_1 = pdf_[number]->xfxQ(id1, x1, Q);
  pdf2_1 = pdf_[number]->xfxQ(id2, x2, Q);

  double weight_0 = pdf1_0 * pdf2_0;
  double weight_1 = pdf1_1 * pdf2_1;
  if ( weight_0 > 0 )
    return weight_1 / weight_0;
  else
    return -1;
}
