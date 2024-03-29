#include "Temp_common.h"

std::string random_string( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}



std::string random_string_restricted( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char {
    const char charset[] = "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}


std::string GAP_logical(bool const& x)
{
  if (x)
    return "true";
  return "false";
}

std::vector<int> StdVectorFirstNentries(size_t const& N)
{
  std::vector<int> eList(N);
  for (size_t i=0; i<N; i++)
    eList[i]=i;
  return eList;
}

void WriteVectorInt_GAP(std::ostream &os, std::vector<int> const& OneInc)
{
  size_t siz=OneInc.size();
  os << "[";
  for (size_t i=0; i<siz; i++) {
    if (i>0)
      os << ",";
    size_t eVal=OneInc[i]+1;
    os << eVal;
  }
  os << "]";
}

std::vector<int> DivideListPosition(int const& len, int const& nbBlock)
{
  std::vector<int> ListVal;
  for (int i=0; i<=nbBlock; i++) {
    double pos_d = (double(i) / double(nbBlock)) * double(len);
    int pos_i;
    NearestInteger_double_int(pos_d, pos_i);
    pos_i = std::max(0, pos_i);
    pos_i = std::min(len, pos_i);
    ListVal.push_back(pos_i);
  }
  return ListVal;
}


