/*
 * FileMgr.cpp
 *
 *  Created on: Jan 15, 2010
 *      Author: mlandis
 */

#include "FileMgr.h"

FileMgr::FileMgr(void) {
	setFileName("");
	setFilePath("");
	setCurrDir(findCurrDir());
	setFilePath(getCurrDir());
}

FileMgr::FileMgr(std::string s) {
	setFileName("");
	setFilePath("");
	parsePathFileNames(s);
	setCurrDir(findCurrDir());
	if (getFilePath() == "") {
		setFilePath(getCurrDir());
	}
}

FileMgr::~FileMgr(void) {
}

std::string FileMgr::getFilePathName(void) {
	std::string delimiter = "/";
	return filePath + delimiter + fileName;
}

void FileMgr::closeFile(std::ifstream &strm) {
	strm.close();
}

bool FileMgr::openFile(std::ifstream &strm) {
	std::string delimiter = "/";
	std::string fileNamePath = filePath + delimiter + fileName;
	strm.open(fileNamePath.c_str(), std::ios::in);
	if (!strm) {
		return false;
	}
	return true;
}

bool FileMgr::openFile(std::ofstream &strm) {
	std::string delimiter = "/";
	std::string fileNamePath = filePath + delimiter + fileName;
	strm.open(fileNamePath.c_str(), std::ios::out);
	if (!strm) {
		return false;
	}
	return true;
}

bool FileMgr::testDir(void) {
	return isDirPresent(filePath);
}

bool FileMgr::testFile(void) {
	return isFilePresent(filePath, fileName);
}

bool FileMgr::listDirContents(void) {
	/* open the directory */
	DIR *dir = opendir(filePath.c_str());
	if (!dir) {
		std::cerr << "Could not find path to directory" << std::endl;
		return false;
	}

	/* read the directory's contents */
	struct dirent *dirEntry;
	while ((dirEntry = readdir(dir)) != NULL) {
		std::cout << dirEntry->d_name << std::endl;
	}

	/* close the directory */
	if ( closedir(dir) == -1 ) {
		std::cerr << "Problem closing directory" << std::endl;
		return false;
	}

	return true;
}

bool FileMgr::parsePathFileNames(std::string s) {
	std::string delimiter = "/";

	if (s.length() == 0) {
		// filePath is empty
		return false;
	}

	// place to separate filepath from filename
	int location = s.find_last_of(delimiter);

	if (location == -1) {
		// filePath contains no "/", contains no directories
		fileName = s;
		filePath = "";
	}

	else if (location == (int)s.length() - 1) {
		// filePath contains "/" as last character, contains no fileName
		s.erase(location);
		fileName = "";
		filePath = s;
		std::cerr << "filePath provided but fileName needed." << std::endl;
		return false;
	}

	else {
		// fileName and filePath exist
		fileName = s.substr(location + 1, s.length() - location - 1);
		filePath = s.substr(0, location);
		return true;
	}

	// return false as last resort
	return false;
}

std::string FileMgr::findCurrDir(void) {
	std::string delimiter = "/";
	char cwd[MAX_DIR_PATH + 1];
	if (!getcwd(cwd, MAX_DIR_PATH + 1)) {
		std::cerr << "Problem finding current working directory!" << std::endl;
		return "";
	}

	std::string curdir = cwd;
	if (curdir.at(curdir.length() - 1) == delimiter[0]) {
		// remove final "/" as needed
		curdir.erase(curdir.length() - 1);
	}

	return curdir;
}

bool FileMgr::isDirPresent(const std::string mp) {
		/* attempt to open the directory */
		DIR *dir = opendir(mp.c_str());
		if (!dir) {
			return false;
		}

		/* close the directory */
		if (closedir(dir) == -1) {
			std::cerr << "Problem closing directory" << std::endl;
		}

		return true;
}

bool FileMgr::isFilePresent(const std::string mp, const std::string mf) {
	/* open the directory */
	DIR *dir = opendir(mp.c_str());
	if (!dir) {
		std::cerr << "Could not find path to directory" << std::endl;
		return false;
	}

	/* read the directory's contents */
	struct dirent *dirEntry;
	bool foundFile = false;
	while ((dirEntry = readdir(dir)) != NULL) {
		std::string temp = dirEntry->d_name;
		if (temp == mf) {
			foundFile = true;
		}
	}

	/* close the directory */
	if (closedir(dir) == -1) {
		std::cerr << "Problem closing directory" << std::endl;
		return false;
	}

	return foundFile;

}

std::string FileMgr::readLineFromFile(int lineNum) {

	int lineCount = 1;
	std::string linestring = "";
	std::ifstream seqStream;

	if (openFile(seqStream) == false) {
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
	}

	if (lineNum < 1) std::cout << "Invalid lineNum value = " << lineNum << std::endl;

	// read line
	while (getline(seqStream, linestring).good()) {
		if (lineCount == lineNum) return linestring;
		lineCount++;
	}

	std::cerr << "Could not find lineNum = " << lineNum << std::endl;
	return "";

}
