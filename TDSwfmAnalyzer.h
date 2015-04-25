#include "TApplication.h"
#include "TSystem.h"
#include "TFile.h"
#include "TBrowser.h"
#include <vector>
#include <iostream>

class RBrowser : public TBrowser
{
  
public:
  
  RBrowser();
  ~RBrowser();
};

RBrowser::RBrowser() : TBrowser()
{
}

RBrowser::~RBrowser()
{
  // exit on close window
  gApplication->Terminate();
}
