#ifndef BASIC_STRING_INCLUDE
#define BASIC_STRING_INCLUDE

#include "Temp_common.h"

std::string STRING_GETENV(std::string const& eStr);
bool STRING_IsStringReduceToSpace(std::string const& eStr);
int STRING_GetCharPositionInString(std::string const& eStr, std::string const& eChar);
std::string StringSubstitution(std::string const& FileIN, std::vector<std::pair<std::string,std::string>> const& ListSubst);
bool IsFullyNumeric(std::string const& eStr);
std::string DoubleTo4dot2f(double const& x);
std::string DoubleTo4dot1f(double const& x);
std::string DoubleToString(double const& x);
std::string IntToString(int const & x);
int StringToInt(std::string const& str);
std::string LongToString(long const & x);
int GetNumberDigit(int const& eVal);
std::string StringNumber(int const& nb, int const& nbDigit);
std::string UpperCaseToLowerCase(std::string const& dataIn);
std::string STRING_RemoveSpacesBeginningEnd(std::string const& eStr);

// Example of use eStrA="A B C D" and eStrB=" "
std::vector<std::string> STRING_Split(std::string const& eStrA, std::string const& eStrB);

// The String is supposed to be "str0" + hs0 + "str1" + hs1 + "hs2"
// and we return a standard vector of [hs0, hs1]
std::vector<std::string> STRING_ParseSingleLine(std::string const& strin, std::vector<std::string> const& LStr);

std::vector<int> STRING_Split_Int(std::string const& eStrA, std::string const& eStrB);

// Example of use eStrA="A B C D" and eStrB=" "
// Difference with the above is that splitting ";;" gives you a list
// of entries as {"", "", ""}, i.e. last and first entry and entry in the middle
std::vector<std::string> STRING_Split_Strict(std::string const& eStrA, std::string const& eStrB);

std::string STRING_Replace(std::string const&eStrA, std::string const& eStrB, std::string const& eStrC);

std::vector<std::string> STRING_SplitCharNb(std::string const& str);

std::string FILE_GetExtension(std::string const& eFile);

#endif
