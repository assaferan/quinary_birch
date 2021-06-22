#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ParseNipp.h"

NippEntry
parseNextGenus(std::ifstream& nippFile, const std::string & line)
{
  NippEntry entry;

  std::istringstream line_str(line);
  std::string desc;
  char next_char;

  line_str >> desc;
  assert( desc == "D=");
  line_str >> entry.disc;

  line_str >> next_char;
  assert( next_char == ';');

  line_str >> desc;
  assert( desc == "GENUS#");
  line_str >> entry.genus;
  
  line_str >> next_char;
  assert(next_char == ';');

  line_str >> desc;
  assert(desc == "MASS=");
  line_str >> entry.mass[0];

  line_str >> next_char;
  assert(next_char == '/');
  line_str >> entry.mass[1];

  line_str >> next_char;
  assert(next_char == ';');

  line_str >> desc;
  assert(desc == "HASSE");
  line_str >> desc;
  assert(desc == "SYMBOLS");
  line_str >> desc;
  assert( desc == "ARE");

  short int symb;
  
  while (line_str)    {
    line_str >> symb;
    entry.HasseSymb.push_back(symb);
  }
  
  next_char = nippFile.peek();

  // There might be trouble at the end of the file - !! TODO
  while (next_char != 'D') {
    LatticeRecord lattice;
    for (size_t i = 0; i < LatticeRecord::VecSize; i++)
      nippFile >> lattice.form[i]; 
    next_char = nippFile.get();
    assert(next_char == ';');
    nippFile >> lattice.numAut;
    entry.lattices.push_back(lattice);
    if (nippFile.eof())
      return entry;
    next_char = nippFile.get();
    assert(next_char == '\n');
    next_char = nippFile.peek();
  }
  
  return entry;
}

std::vector<NippEntry>
ParseNipp::parseDisc(const std::string & fname, const Z & disc)
{
  std::vector<NippEntry> genera;
  std::ifstream nippFile(fname);
  std::ostringstream disc_str, find_str;
  
  find_str << "D=";
  disc_str << disc;
  for (size_t i = 0; i < 5 - disc_str.str().size(); i++)
    find_str << " ";
  find_str << disc_str.str();

  if (nippFile.is_open())
    {
      std::string line;
      bool found = false;
      bool done = false;
      while ((!done) && (std::getline(nippFile, line)))
	{
	  size_t start = line.find(find_str.str(), 0);
	  if (start != std::string::npos)
	    {
	      // found it
	      found = true;
	      NippEntry latGen = parseNextGenus(nippFile, line);
	      genera.push_back(latGen);
	    }
	  else
	    if (found) done = true; 
	}
    }
  else throw std::runtime_error("Unable to open nipp file.");
  nippFile.close();
  return genera;
}
