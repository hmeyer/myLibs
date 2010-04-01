
#ifndef TEMPLATE_TO_STRING
#define TEMPLATE_TO_STRING

#include <sstream>

template <class TYPE> std::string toString( const TYPE & t ) {
 std::ostringstream os;
 os << t;
 return os.str();
}

template<typename T>
T fromString(const std::string& s) {
 T toT;
 std::istringstream(s) >> toT;
 return toT;
}

#endif
