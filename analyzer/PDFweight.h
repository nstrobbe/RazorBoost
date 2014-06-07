#ifndef PDFWEIGHT_H
#define PDFWEIGHT_H
//-----------------------------------------------------------------------------
// File:        PDFweight.h
// Created:     05-June-2014 - reworking of some old code
// Author:      Harrison B. Prosper      
//-----------------------------------------------------------------------------
#include <string>
#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"
//-----------------------------------------------------------------------------
using namespace std;
//-----------------------------------------------------------------------------
class PDFweight
{
 public:
  PDFweight() {}
  PDFweight(std::string pdfName, std::string pdfNameOrig="");
  ~PDFweight();

  double operator()(int id1,     /// PDGID of parton 1
		    int id2,     /// PGDID of parton 2
		    double x1, double x2, double Q, 
		    int number); /// PDF member number

 private:
  LHAPDF::PDFSet pdfset_;
  std::vector<LHAPDF::PDF*> pdf_;
  LHAPDF::PDF* pdf0_;
};

#endif

