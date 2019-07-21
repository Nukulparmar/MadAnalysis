#include "SampleAnalyzer/User/Analyzer/cms_sus_18_002.h"
#include "SampleAnalyzer/Process/Analyzer/AnalyzerManager.h"
#include "SampleAnalyzer/Commons/Service/LogStream.h"

// -----------------------------------------------------------------------------
// BuildTable
// -----------------------------------------------------------------------------
void BuildUserTable(MA5::AnalyzerManager& manager)
{
  using namespace MA5;
  manager.Add("cms_sus_18_002",new cms_sus_18_002);
}
