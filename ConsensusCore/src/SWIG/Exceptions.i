
%{
/* Includes the header in the wrapper code */
#include <iostream>
#include <ConsensusCore/Types.hpp>
using namespace ConsensusCore;
%}

%include <ConsensusCore/Types.hpp>

%include <exception.i>

%exception {
  try {
    $action
  } catch (const std::exception& e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const ConsensusCore::ErrorBase& e) {
      SWIG_exception(SWIG_RuntimeError, e.Message().c_str());
  } catch (const ConsensusCore::ExceptionBase& e) {
      std::cout << "[Undeclared exception thrown---amend C++ code to declare]" << std::endl;
      SWIG_exception(SWIG_RuntimeError, e.Message().c_str());
  }
}
