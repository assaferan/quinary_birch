#ifndef BASIC_FILE_INCLUDE
#define BASIC_FILE_INCLUDE

#include "Temp_common.h"
#include "Basic_string.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


void CopyOperation(std::string const& SrcFile, std::string const& DstFile);
bool IsExistingFile(std::string const& eFile);
void IsExistingFileDie(std::string const& eFile);
std::vector<std::string> ReadFullFile(std::string const& eFile);
std::string FILE_GetNakedFilename(std::string const& eFileFull);
std::string FILE_RemoveExtension(std::string const& eFileIn);
std::string FILE_GetDirectoryOfFileName(std::string const& eFileFull);
std::vector<std::string> FILE_GetDirectoryListFile(std::string const& eDir);
bool FILE_IsDirectoryEmpty(std::string const& eDir);
bool FILE_CheckFinalShashDirectory(std::string const& eDir);
bool FILE_IsRegularFile(std::string const& eFile);
std::vector<std::string> FILE_GetDirectoryFilesRecursively(std::string const& eDir);
std::vector<std::string> FILE_DirectoryFilesSpecificExtension(std::string const& ePrefix, std::string const& eExtension);
std::vector<std::string> FILE_DirectoryMatchingPrefixExtension(std::string const& ePrefix, std::string const& eExtension);

// If the ePrefix ends with / then we do recursive listing
// If the output is of the form path/WWM_output_
// then returns all the files having the 
std::vector<std::string> FILE_DirectoryFilesSpecificExtension_Gen(std::string const& ePrefix, std::string const& eExtension);
bool IsExistingDirectory(std::string const& ThePrefix);
void RemoveEmptyDirectory(std::string const& eDir);
void RemoveFile(std::string const& eFile);
void RemoveFileIfExist(std::string const& eFile);
std::string FILE_RemoveEndingExtension(std::string const& FileName, std::string const& TheExtension);
void RemoveFileSpecificExtension(std::string const& ThePrefix, std::string const& TheExtension);
void RemoveFileInDirectory(std::string const& ThePrefix);
int FILE_GetNumberLine(std::string const& eFile);

#ifndef WINDOWS
std::string GetCurrentDirectory();
#endif

#ifndef WINDOWS
std::string FILE_GetAbsoluteDirectory(std::string const& ePrefix);
#endif

bool FILE_CheckPrefix(std::string const& ePrefix);
std::string ExtractDirectoryFromFileString(std::string const& eFile);
bool FILE_IsFileMakeable(std::string const& eFile);
void CreateDirectory(std::string const& eDir);
void CreateDirectory_V1(std::string const& eDir);
std::vector<std::string> ls_operation(std::string const& ThePrefix);

struct TempFile {
private:
  std::string FileName;
public:
  TempFile() = delete;
  TempFile(std::string const& eFile);
  TempFile(char* eFile);
  ~TempFile();
  //
  bool IsExisting() const; 
  std::string string() const;
};

// This does not provide all the facility that we are after
// The problem is that when we have an exit(1)
// the destructors of the variables are not called.
//
// This means that we need to put a throw statement.
// in order to have the temporary directory eliminated.
// So eliminate all the exit(1) statement and replace by throw
// and catch the throws.
struct TempDirectory {
private:
  std::string DirName;
  bool IsInitialized;
public:
  TempDirectory();
  TempDirectory(std::string const& eDir);
  TempDirectory& operator=(TempDirectory&& eTemp);
  TempDirectory(TempDirectory && eTemp);
  TempDirectory(const TempDirectory & eTemp) = delete;
  TempDirectory& operator=(const TempDirectory & eTemp) = delete;
  
  ~TempDirectory();
  //
  bool IsExisting() const;
  std::string str() const;
};


struct CondTempDirectory {
private:
  bool used;
  std::string DirName;
public:
  CondTempDirectory();
  CondTempDirectory(bool const& eUsed, std::string const& eDir);
  CondTempDirectory& operator=(CondTempDirectory&& eTemp);
  CondTempDirectory(CondTempDirectory && eTemp);
  ~CondTempDirectory();
  //
  bool IsExisting() const;
  bool usedness() const;
  std::string str() const;
};



#endif
