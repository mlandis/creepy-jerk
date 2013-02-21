/*
 * FileMgr.h
 *
 *  Created on: Jan 15, 2010
 *      Author: mlandis
 */

#ifndef FILEMGR_H_
#define FILEMGR_H_
#define MAX_DIR_PATH 2048

#include <dirent.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <unistd.h>

class FileMgr {
public:
				FileMgr(void);
				FileMgr(std::string s);
	virtual 	~FileMgr(void);
	std::string	getCurrDir(void)				{ return currDir; }
	std::string	getFileName(void)				{ return fileName; }
	std::string	getFilePath(void)				{ return filePath; }
	std::string	getFilePathName(void);
	std::string	readLineFromFile(int);
	void		setCurrDir(std::string s)			{ currDir = s; }
	void		setFileName(std::string s)			{ fileName = s; }
	void		setFilePath(std::string s)			{ filePath = s; }
	void		closeFile(std::ifstream &strm);
	bool		openFile(std::ifstream &strm);
	bool		openFile(std::ofstream &strm);
	bool		testDir(void);
	bool		testFile(void);
	bool		listDirContents(void );
	bool		parsePathFileNames(std::string s);

private:
	std::string	findCurrDir(void);
	bool		isDirPresent(const std::string mp);
	bool		isFilePresent(const std::string mp, const std::string mf);
	std::string	currDir;
	std::string	filePath;
	std::string	fileName;

};

#endif /* FILEMGR_H_ */
