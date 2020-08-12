import logging
from pysnptools.util.mapreduce1.runner import *
import os
import cPickle as pickle
import subprocess, sys, os.path
import multiprocessing
import pysnptools.util as pstutil
import pdb
from collections import defaultdict
import tarfile as tarfileLibrary
import ctypes
import datetime
from tempfile import TemporaryFile

class Hadoop(Runner):

    '''
    Old code to run on Hadoop. Not currently supported.
    '''


    fileshare = "/user"

    def __init__(self, taskcount, mapmemory=-1, reducememory=-1, mkl_num_threads=None,queue="default",skipsourcecheck=False,skipdatacheck=False,logging_handler=logging.StreamHandler(sys.stdout)):
        logger = logging.getLogger()
        if not logger.handlers:
            logger.setLevel(logging.INFO)
        for h in list(logger.handlers):
            logger.removeHandler(h)
        logger.addHandler(logging_handler)
        if logger.level == logging.NOTSET:
            logger.setLevel(logging.INFO)

        logging.info('constructing Hadoop runner')
        self.taskcount = taskcount
        self.mapmemory = mapmemory
        self.reducememory = reducememory
        self.mkl_num_threads = mkl_num_threads
        self.queue = queue
        self.skipsourcecheck = skipsourcecheck
        self.skipdatacheck = skipdatacheck
        logging.info('done constructing Hadoop runner')

    def run(self, distributable):
        logging.info('Hadoop runner is running a distributable')

        # Check that the local machine has python path set
        localpythonpath = os.environ.get("PYTHONPATH") #!!should it be able to work without pythonpath being set (e.g. if there was just one file)? Also, is None really the return or is it an exception.
        if localpythonpath == None: raise Exception("Expect local machine to have 'pythonpath' set")

        remotewd, run_dir_abs, run_dir_rel = self.create_run_dir()

        result_local = os.path.join(run_dir_rel,"result.p")
        result_hdfs = "hdfs:" + os.path.join(run_dir_abs,"result.p").replace("\\","/")
        result_remote = os.path.join(distributable.tempdirectory,"result.p")

        fileInWorkingDirectoryList =[]

        tgzListPythonSettings = self.FindOrCreatePythonSettings(remotewd)

        remotepythonpath, tgzListPythonPath = self.FindOrCreateRemotePythonPath(localpythonpath, remotewd)

        tgzList = []
        inputOutputCopier = HadoopCopier(remotewd, fileInWorkingDirectoryList, tgzList, self.skipdatacheck) #Create the object that copies input and output files to where they are needed
        inputOutputCopier.input(distributable) # copy of the input files to where they are needed (i.e. the cluster)


        batfilename_abs_list = self.create_bat_file(distributable, remotepythonpath, remotewd, run_dir_abs, run_dir_rel, result_remote, result_hdfs)

        self.submit_to_cluster(batfilename_abs_list, fileInWorkingDirectoryList, tgzList, tgzListPythonPath, tgzListPythonSettings, distributable, remotewd, run_dir_abs, run_dir_rel)

        inputOutputCopier.output(distributable) # copy the output file from where they were created (i.e. the cluster) to the local computer

        subprocess.check_output("%HADOOP_HOME%\\bin\\hadoop fs -copyToLocal {0} {1}\n".format(result_hdfs, result_local),stderr=subprocess.STDOUT,shell=True)
        with open(result_local, mode='rb') as f:
            result = pickle.load(f)

        #logging.info('Done: Hadoop runner is running a distributable. Returns {0}'.format(result))
        return result

    def submit_to_cluster(self, batfilename_rel_list, fileInWorkingDirectoryList, tgzList, tgzListPythonPath, tgzListPythonSettings, distributable, remotewd, run_dir_abs, run_dir_rel):
        logging.info('Hadoop runner is submitting to cluster')

        #!! e.g. hdfs://rr1-n13-02-c02/user/carlk/inputs.tgz#inputs,hdfs://rr1-n13-02-c02/user/carlk/datasets.tgz#datasets,hdfs://rr1-n13-02-c02/user/carlk/src.tgz#src
        #!! could do this functionally
        archivesStringList = []
        for tgz in tgzList:
            archiveString = "hdfs:{0}#{1}".format(tgz[1],os.path.splitext(tgz[0])[0])
            archivesStringList.append(archiveString)
        archivesStringList.append(tgzListPythonSettings)
        for tgz in tgzListPythonPath:
            archivesStringList.append(tgz)

        #e.g. distMapper.bat,distReducer.bat
        filesString = ",".join(batfilename_rel_list+fileInWorkingDirectoryList)

        taskIndexDir = run_dir_rel + os.path.sep + "input"
        pstutil.create_directory_if_necessary(taskIndexDir,isfile=False)

        #zgoal = int(SP.ceil(SP.log(self.taskcount)/SP.log(10)))
        with open(taskIndexDir +  os.path.sep + "taskIndexList.txt","w") as taskIndexListFile:
            for taskIndex in xrange(self.taskcount):
                taskIndexListFile.write("{0}\n".format(taskIndex)) # str(taskIndex).zfill(zgoal)))

        #hadoop fs -rmr runs/2013-08-02_13_51_42
        #hadoop fs -copyFromLocal runs\2013-08-02_13_51_42 runs/2013-08-02_13_51_42
        #hadoop jar %HADOOP_HOME%\lib\hadoop-streaming.jar ^
        #        -archives "hdfs:/user/carlk/source/carlkextranet05312013/ERG01/src/tests/datasets.2013-07-31_11_12_11.tgz#datasets,hdfs:/user/carlk/runs/pythonpath.0.src.2013-07-31_14_30_56/src.tgz#pythonpath.0.src" ^
        #        -files "hdfs:/user/carlk/runs/2013-08-02_13_51_42/distMapper.bat,hdfs:/user/carlk/runs/2013-08-02_13_51_42/distReducer.bat,hdfs:/user/carlk/runs/2013-08-02_13_51_42/distributable.p" ^
        #        -input "runs/2013-08-02_13_51_42/input" ^
        #        -output "runs/2013-08-02_13_51_42/output" ^
        #        -mapper "distMapper.bat" ^
        #       -reducer "distReducer.bat"
        #hadoop fs -cat runs/2013-08-02_13_51_42/output/part-00000  | more
        s00 = r"%HADOOP_HOME%\bin\hadoop fs -rmr -skipTrash {0}".format(run_dir_rel.replace("\\","/"))
        s0 = r"%HADOOP_HOME%\bin\hadoop fs -copyFromLocal {0} {1}".format(run_dir_rel, run_dir_rel.replace("\\","/"))


        #-D mapreduce.reduce.shuffle.connect.timeout=3600000 ^
        #-D io.sort.mb=1400 ^
        #-D job.end.retry.interval=3600000 ^
        #-D mapred.tasktracker.expiry.interval=3600000 ^

        logging.info("running {0}".format(str(distributable)))
        

        s = r"""%HADOOP_HOME%\bin\hadoop jar %HADOOP_HOME%\lib\hadoop-streaming.jar ^
        -archives "{0}" ^
        -files "{1}" ^
        -D mapred.job.name="{8}" ^
        -D mapred.map.tasks={4} ^
        -D mapred.reduce.tasks=1 ^
        -D mapred.job.map.memory.mb={5} ^
        -D mapred.job.reduce.memory.mb={6} ^
        -D mapred.task.timeout={7} ^
        -D mapred.job.queue.name="{9}" ^
        -input {2} ^
        -output {3} ^
        -inputformat org.apache.hadoop.mapred.lib.NLineInputFormat ^
        -mapper "distMapper.bat" ^
        -reducer "distReducer.bat"
            """.format(
                    ",".join(archivesStringList),       #0
                    filesString,                        #1
                    taskIndexDir.replace("\\","/"),     #2
                    (run_dir_rel + os.path.sep + "output").replace("\\","/"), #3
                    self.taskcount,                     #4
                    self.mapmemory,                     #5
                    self.reducememory,                  #6
                    0,                                  #7
                    str(distributable),                 #8
                    self.queue                         #9
                    )
        runHadoopFileName = run_dir_rel + os.path.sep + "runHadoop.bat"
        logging.info("Hadoop runner is creating '{0}'".format(runHadoopFileName))
        with open(runHadoopFileName, "w") as runHadoopFile:
            runHadoopFile.write("call {0}\n".format(s00))
            runHadoopFile.write("call {0}\n".format(s0))
            runHadoopFile.write("call {0}\n".format(s))

        sOneLine = "".join(s.split("^\n"))

        
        logging.info("Hadoop runner running the copyFromLocal")
        with TemporaryFile() as output:
            stdout0 = subprocess.check_output(s0,stderr=output,shell=True)
            output.seek(0)
            stderr0 = output.read()
        logging.info("Result from 'Hadoop runner running the copyFromLocal' is stdout='{0}', stderr='{1}'".format(stdout0, stderr0))
        if stderr0 != "" : raise Exception("Stderr from command: '{0}'".format(stderr0))
        logging.info("Hadoop runner running the streamingjar")

        with TemporaryFile() as output:
            stdout = subprocess.check_output(sOneLine,stderr=output,shell=True)
            output.seek(0)
            stderr = output.read()
        logging.info("Result from 'Hadoop runner running the streamingjar' is stdout='{0}', stderr='{1}'".format(stdout, stderr))
        logging.info('Done: Hadoop runner is submitting to cluster')
        #if stderr != "" : raise Exception("Stderr from command: '{0}'".format(stderr))

    def FindOrCreatePythonSettings(self, remotewd):
        localpythonpathsetting = r"\\msr-arrays\scratch\msr-pool\eScience3\.continuum" # os.path.join(os.environ.get("userprofile"),".continuum")
        lastFolderName = os.path.split(os.path.normpath(localpythonpathsetting))[1]
        #pstutil.create_directory_if_necessary(localpythonpathsetting,isfile=False)

        #Must set assume_changed=True for otherwise hidden .continuum file to be used.
        tgzName = HadoopCopier.CheckUpdateTgz(localpythonpathsetting, subsubItemList1=None, skipcheck=False, filter_hidden=False)
        hdfstgz = "hdfs:{3}/{2}.{1}/{2}.tgz".format(None,str(datetime.datetime.fromtimestamp(os.path.getmtime(tgzName)))[:19].replace(" ","_").replace(":","_"),lastFolderName,remotewd)
        Hadoop.hdfsCopyFromLocalIfNotThere(tgzName, hdfstgz)
        return hdfstgz + "#" + lastFolderName

    def create_distributablep(self, distributable, run_dir_abs, run_dir_rel):
        logging.info('Hadoop runner is pickling distributable')
        distributablep_filename_rel = os.path.join(run_dir_rel, "distributable.p")
        #distributablep_filename_abs = os.path.join(run_dir_abs, "distributable.p")
        pstutil.create_directory_if_necessary(distributablep_filename_rel)
        with open(distributablep_filename_rel, mode='wb') as f:
            pickle.dump(distributable, f, pickle.HIGHEST_PROTOCOL)
        logging.info('Done: Hadoop runner is pickling distributable')
        return distributablep_filename_rel

    def create_bat_file(self, distributable, remotepythonpath, remotewd, run_dir_abs, run_dir_rel, result_remote, result_hdfs):
        logging.info('Hadoop runner is creating bat file')

        outFileList = Hadoop.RecursivelyGetAllOutputs(distributable)

        distributablep_filename_rel = self.create_distributablep(distributable, run_dir_abs, run_dir_rel)

        distributable_py_file = os.path.join(os.path.dirname(__file__),"..","distributable.py")
        if not os.path.exists(distributable_py_file): raise Exception("Expect file at " + distributable_py_file + ", but it doesn't exist.")
        localfilepath, file = os.path.split(distributable_py_file)
        remoteexepath = os.path.join(remotepythonpath.split(';')[0],"fastlmm","util") #!!shouldn't need to assume where the file is in source

        batfilename_abs_list = []
        for part in ["Mapper","Reducer"]:
            command_string = remoteexepath + os.path.sep + file + r""" distributable.p "Local{0}({1},""{2}"",mkl_num_threads={3},logging_handler=logging.StreamHandler())" """.format(
                    part,
                    self.taskcount,
                    result_remote.replace("\\","/"), #change to DOS separator to Unix separator because python will work with either and this avoid problems with parsing the batch file
                    self.mkl_num_threads)
            batfilename_rel = os.path.join(run_dir_rel,"dist{0}.bat".format(part))
            batfilename_abs = "hdfs:" + os.path.join(run_dir_abs,"dist{0}.bat".format(part)).replace("\\","/")
            batfilename_abs_list.append(batfilename_abs)
            pstutil.create_directory_if_necessary(batfilename_rel, isfile=True)
            with open(batfilename_rel, "w") as batfile:
                batfile.write("@set path={0};{0}\Scripts;%path%\n".format(r"c:\GCD\esciencepy"))
                batfile.write("@set PYTHONPATH={0}\n".format(remotepythonpath))
                batfile.write("@set home=%cd%\n")
                #batfile.write("@mklink /d .continuum continuum\n")
                #batfile.write("@dir /s\n")
                #batfile.write("@set R_HOME={0}\n".format(os.path.join(remotepythoninstall,"R-2.15.2")))
                #batfile.write("@set R_USER={0}\n".format("."))
                batfile.write("@mkdir {0}\n@mkdir {0}\\tex.cache\n@set MPLCONFIGDIR={0}\n".format(".matplotlib"))
                batfile.write("@mkdir {0}\nset IPYTHONDIR={0}\n".format(".ipython"))
                #batfile.write("xcopy /d /e /s /c /h /i continuum .continuum\n")
                batfile.write("@call python {0}\n".format(command_string))
                if part == "Reducer":
                    batfile.write("@call %HADOOP_HOME%\\bin\\hadoop fs -rm {0} -skipTrash\n".format(result_hdfs))
                    batfile.write("@call %HADOOP_HOME%\\bin\\hadoop fs -copyFromLocal {0} {1}\n".format(result_remote, result_hdfs))
                    for outfile in outFileList:
                        hdfsOutFile = remotewd + "/" + outfile
                        batfile.write("@call %HADOOP_HOME%\\bin\\hadoop fs -rm {0}\n".format(hdfsOutFile))
                        batfile.write("@call %HADOOP_HOME%\\bin\\hadoop fs -copyFromLocal {0} {1}\n".format(outfile, hdfsOutFile))
        picklefilename_abs = "hdfs:" + os.path.join(run_dir_abs,"distributable.p").replace("\\","/")
        batfilename_abs_list.append(picklefilename_abs)
        logging.info('Done: Hadoop runner is creating bat file')
        return batfilename_abs_list

    @staticmethod
    def RecursivelyGetAllOutputs(item):
        outputList = []
        ListCopier([],outputList).output(item)
        return outputList

    def FindOrCreateRemotePythonPath(self, localpythonpath, remotewd):
        #input: something like: 'D:\\Source\\carlkextranet05312013\\ERG01\\src\\' and maybe a second item, etc
        #sideeffect: create D:\\Source\\carlkextranet05312013\\ERG01\\src.tgz and 2nd item
        #sideeffect: copy if newer to to hdfs /user/carlk/runs/pythonpath.src.0.<moddate>/src.tgz and ....1.....2nditem
        #return 1 list of "hdfs /user/carlk/runs/pythonpath.0.<moddate>/src.tgz#pythonpath.0.src"
        #       2 remotepythonpath, e.g. "pythonpath.0.src;pythonpath.1.2nditem"

        remotepythonpath_list = []
        tgzList = []
        for i, localpythonpathdir in enumerate(localpythonpath.split(';')):
            tgzName = HadoopCopier.CheckUpdateTgz(localpythonpathdir,skipcheck=self.skipsourcecheck)
            lastFolderName = os.path.split(os.path.normpath(localpythonpathdir))[1]
            hdfstgz = "hdfs:{3}/pythonpath.{0}.{2}.{1}/{2}.tgz".format(i,str(datetime.datetime.fromtimestamp(os.path.getmtime(tgzName)))[:19].replace(" ","_").replace(":","_"),lastFolderName,remotewd)
            Hadoop.hdfsCopyFromLocalIfNotThere(tgzName, hdfstgz)
            remotepythonpathdir = "pythonpath.{0}.{1}".format(i,lastFolderName)
            remotepythonpath_list.append(remotepythonpathdir)
            tgzList.append(hdfstgz + "#" + remotepythonpathdir)

        remotepythonpath = ";".join(remotepythonpath_list)
        
        #these are parallel
        return remotepythonpath, tgzList

    @staticmethod
    def hdfsCopyFromLocalIfNotThere(tgzName, hdfstgz):
        # if it is there won't copy. If it isn't there will copy.
        subprocess.check_output(r"type {0} | %HADOOP_HOME%\bin\Hadoop fs -put - {1}".format(tgzName, hdfstgz[5:]),stderr=subprocess.STDOUT,shell=True)

    def create_run_dir(self):
        username = os.environ["USERNAME"]
        localwd = os.getcwd()
        if localwd.startswith("\\\\"):
            remotewd = self.fileshare + os.path.sep + username + os.path.sep + "\\".join(localwd.split('\\')[4:])
        else:
            remotewd = self.fileshare + os.path.sep + username + os.path.splitdrive(localwd)[1]  #using '+' because 'os.path.join' isn't work with shares
        remotewd = remotewd.replace("\\","/")
        if remotewd.endswith("/"): # remove trailing /
            remotewd = remotewd[:-1]
        run_dir_rel = os.path.join("runs",pstutil._datestamp(appendrandom=True))
        pstutil.create_directory_if_necessary("runs",isfile=False)
        if not os.path.isfile(".ignoreTgzChange"):
            with open("runs" +  os.path.sep + ".ignoreTgzChange","w") as ignoreFile:
                ignoreFile.write("\n")


        run_dir_abs = "/user/{0}/{1}".format(username,run_dir_rel)
        #!! hadoop_create_directory_if_necessary(run_dir_abs,isfile=False)
        return remotewd, run_dir_abs, run_dir_rel

    #!! move these hadoop commands to a library
    @staticmethod
    def hadoop_create_directory_if_necessary(name, isfile=True):    
        import os
        if isfile:
            directory_name = os.path.dirname(name)
        else:
            directory_name = name
        hadoop_makedirs(directory_name)

    #!! what if already is there?
    @staticmethod
    def hadoop_makedirs(directory_name):
        hadoop_command("fs -mkdir {0}".format(directory_name))

    @staticmethod
    def hadoop_command(command_string):
        rc = os.system(os.environ.get("HADOOP_HOME") + os.path.sep + "bin" + os.path.sep + "hadoop " + command_string)
        if rc != 0: raise Exception("Hadoop command '{0}' fails with rc={1}.".format(command_string,rc))


class ListCopier(object): #Implements ICopier

    def __init__(self, inputList, outputList): #The list will be modified
        if len(inputList) != 0 : raise Exception("Expect inputList to start empty")
        if len(outputList) != 0 : raise Exception("Expect outputList to start empty")
        self.inputList = inputList
        self.outputList = outputList

    def input(self,item):
        if isinstance(item, str):
            self.inputList.append(item)
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        #else do nothing
            pass # ignore

    def output(self,item):
        if isinstance(item, str):
            self.outputList.append(item)
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        #else do nothing
            pass # ignore

class HadoopCopier(object): #Implements ICopier

    def __init__(self, remotewd, fileInWorkingDirectoryList, tgzList, skipdatacheck): #The two lists will be modified
        if len(fileInWorkingDirectoryList) != 0 : raise Exception("Expect fileInWorkingDirectoryList to start empty")
        if len(tgzList) != 0 : raise Exception("Expect tgzList to start empty")
        self.remotewd = remotewd
        self.fileInWorkingDirectoryList = fileInWorkingDirectoryList
        self.tgzList = tgzList
        self.skipdatacheck = skipdatacheck


    def input(self,item):
        inputList = self.RecursivelyGetAllInputs(item)

        #group by subfolder, treating files in the working directory as special and giving an error with higher levels
        fileInWorkingDirectoryList, subDirectoryToSubSubItemList = self.GroupByTopLevelSubFolder(inputList)

        #create or update a tgz for each directory
        for directory,subsubItemList in subDirectoryToSubSubItemList.iteritems():
            tgzName = HadoopCopier.CheckUpdateTgz(directory, subsubItemList1=subsubItemList, skipcheck=self.skipdatacheck)
            hdfsName = "{0}/{1}.{2}.tgz".format(self.remotewd,directory,str(datetime.datetime.fromtimestamp(os.path.getmtime(tgzName)))[:19].replace(" ","_").replace(":","_")).replace("\\","/")
            self.tgzList.append((tgzName,hdfsName))
            hdfsNameWild = "{0}/{1}.{2}.tgz".format(self.remotewd,directory,'*').replace("\\","/")
            lookout = subprocess.check_output(r"%HADOOP_HOME%\bin\Hadoop fs -ls {0}".format(hdfsName),stderr=subprocess.STDOUT,shell=True)
            if "No such file or directory" in lookout:
                subprocess.check_output(r"%HADOOP_HOME%\bin\Hadoop fs -rm {0}".format(hdfsNameWild),stderr=subprocess.STDOUT,shell=True)
                subprocess.check_output(r"%HADOOP_HOME%\bin\Hadoop fs -copyFromLocal {0} {1}".format(tgzName, hdfsName),stderr=subprocess.STDOUT,shell=True)

        for file in fileInWorkingDirectoryList:
            if os.path.getsize(file) > 10 * 1024 * 1024:
                logging.warn("File '{0}' is large ({1}GB) and would be more efficiently  distributed on Hadoop if placed in subfolder of the working directory.".format(file, round(os.path.getsize(file)/(1024*1024))))
            self.fileInWorkingDirectoryList.append(file)
            # filesInWorkingDirectory copy to hdfs

    @staticmethod
    def CheckUpdateTgz(directory0, subsubItemList1 = None, skipcheck=False, filter_hidden=True):
        if subsubItemList1 == None:
                subsubItemList1 = [[]]
        directory = os.path.normpath(directory0)
        tgzName =  directory + ".tgz"

        if not skipcheck and os.path.exists(tgzName) and not [] in subsubItemList1 and not filter_hidden:
            if not os.path.isfile(tgzName): raise Exception("Expect '{0}' to be a file.".format(tgzName))
            logging.info("Making list of any files already in {0}".format(tgzName))
            tarfile = tarfileLibrary.open(tgzName, "r")
            for tarinfo in tarfile.getmembers():
                if tarinfo.isfile():
                    filenamein = tarinfo.filename.replace("/",os.path.sep)  #"/" is the tar separator
                    filename = os.path.join(directory,filenamein)
                    if os.path.exists(filename) and os.path.isfile(filename):
                        subsubItemList1.append(filenamein.split("/"))  #"/" is the tar separator
            tarfile.close()

        subsubItemListOut = HadoopCopier.RemoveRedundant(subsubItemList1)

        if not HadoopCopier.IsTarFileUpToDate(directory, subsubItemListOut, tgzName, filter_hidden):
            HadoopCopier.CreateNewTarFile(directory, subsubItemListOut, tgzName, filter_hidden)
        return tgzName

    @staticmethod
    def CreateNewTarFile(directory, subsubItemList, tgzName, filter_hidden=True):

        # logging.info("{0}, {1}, {2}".format(directory, subsubItemList, tgzName))
        directory1 = os.path.normpath(directory)

        tgzNameTemp = tgzName + ".tar" #create a temp file so that both files will exist for a while
        logging.info("Creating '{0}'".format(tgzNameTemp))
        tarfile = tarfileLibrary.open(tgzNameTemp,"w:")
        for subsubItem in subsubItemList:
            tarName =  "/".join(subsubItem)  # tar files use "/" not os.path.sep
            tarfile.add(os.path.normpath(directory + "/" + tarName), tarName, recursive=True, filter=lambda x, directory1=directory1,filter_hidden=filter_hidden : HadoopCopier.tarfile_filter_hidden(x,directory1,filter_hidden))
        tarfile.close()
        logging.info("Compressing '{0}'".format(tgzNameTemp))
        subprocess.call(r"c:\cygwin64\bin\gzip.exe --force --fast {0}".format(tgzNameTemp), shell=True)
        logging.info("Finished Compressing '{0}'".format(tgzNameTemp))
        if os.path.exists(tgzName):
            os.remove(tgzName)
        os.rename(tgzNameTemp+".gz",tgzName)


    @staticmethod
    def RemoveRedundant(subsubItemList1):
        subsubItemListOut = []
        for index1,subsubItem1 in enumerate(subsubItemList1):
            if HadoopCopier.shouldinclude(index1,subsubItem1,subsubItemList1):
                subsubItemListOut.append(subsubItem1)
        return subsubItemListOut

    @staticmethod
    def shouldinclude(index1,subsubItem1,subsubItemList1):
        for index2,subsubItem2 in enumerate(subsubItemList1):
            if HadoopCopier.knocksout(index1,subsubItem1,index2,subsubItem2):
                return False
        return True

    @staticmethod
    def knocksout(index1,subsubItem1,index2,subsubItem2):
        if index1 == index2:
            return False
        if len(subsubItem1) == len(subsubItem2) and index1 > index2:
            return False
        if len(subsubItem1) >= len(subsubItem2):
            for i in xrange(len(subsubItem2)):
                if subsubItem1[i] != subsubItem2[i]:
                    return False
            return True
        return False

    def RecursivelyGetAllInputs(self, item):
        inputList = []
        ListCopier(inputList,[]).input(item)
        return inputList

    def GroupByTopLevelSubFolder(self, inputList):
        filesInWorkingDirectory = []
        subDirectoryToSubSubItemList = defaultdict(list)
        for item in inputList:
            normpath = os.path.normpath(item)
            try:
                relpath = os.path.relpath(normpath)
            except :
                raise Exception("for Hadoop input files must be in or below the current directory ('{0}' is not)".format(item))
            
            if not os.path.exists(relpath): raise Exception("Expect input file (or directory) '{0}' but it doesn't exist".format(relpath))
            parts = relpath.split(os.path.sep)
            if os.path.isfile(relpath): #If the input item is a file ...
                if len(parts) == 1: #If the input item is in the working directory
                    filesInWorkingDirectory.append(relpath)  #add it to the filesInWorkingDirectory list
                else: #A file in a subfolder
                    subDirectoryToSubSubItemList[parts[0]].append(parts[1:]) #Add it to the list for the subfolder
            else: #A folder
                if not os.path.isdir(relpath): raise Exception("assert")
                if len(parts) == 1:
                    subDirectoryToSubSubItemList[relpath].append([])
                else:
                    subDirectoryToSubSubItemList[parts[0]].append(parts[1:])
        return filesInWorkingDirectory, subDirectoryToSubSubItemList
    
    @staticmethod
    def IsTarFileUpToDate(directory, subsubItemList, tgzName, filter_hidden):
        if os.path.exists(tgzName):
            logging.info("Checking if '{0}' is up to date".format(tgzName))
            if not os.path.isfile(tgzName): raise Exception("Expect '{0}' to be a file.".format(tgzName))
            tarfile = tarfileLibrary.open(tgzName, "r")
            isUpToDate = HadoopCopier.IsTarFileUpToDateInternal(directory, subsubItemList, tarfile, filter_hidden)
            tarfile.close()
            logging.info("'{0}' up to date? {1} ".format(tgzName, isUpToDate))
            return isUpToDate
        else:
            return False


    @staticmethod
    def IsTarFileUpToDateInternal(directory, subsubItemList, tgzFile, filter_hidden):
        howToIgnoreString = "To ignore changes in a directory add a file '{0}'".format(HadoopCopier._ignoreTgzChangeFileName)

        for subsubItem in subsubItemList:
            tarName =  "/".join(subsubItem) # use "/" instead of os.path.sep because that is what the tar file uses
            winfileOrDirectory = os.path.normpath(directory + os.path.sep + tarName)
            try:
                member = tgzFile.getmember(tarName)
            except Exception as e:
                logging.info("'{0}' not up to date because of exception {1}. ({2})".format(tarName, e, howToIgnoreString))
                return False;
            else:
                if os.path.isfile(winfileOrDirectory):
                    if int(os.path.getmtime(winfileOrDirectory)) != int(member.mtime):
                        logging.info("'{0}' not up to date because of date change in '{1}'. ({2})".format(tgzFile, winfileOrDirectory, howToIgnoreString))
                        return False;
                else:
                    if not os.path.isdir(winfileOrDirectory): raise Exception("Expect '{0}' to be a directory. ({1})".format(winfileOrDirectory, howToIgnoreString))
                    for winfileOrDirectory2 in HadoopCopier.WalkWithoutHidden(winfileOrDirectory, filter_hidden):
                        tarName2 = os.path.relpath(winfileOrDirectory2, winfileOrDirectory).replace('\\','/')
                        try:
                            member = tgzFile.getmember(tarName2)
                        except KeyError:
                            logging.info("'{0}' not up to date because '{1}' is not found. ({2})".format(tgzFile, tarName2, howToIgnoreString))
                            return False;
                        else:
                            if os.path.isfile(winfileOrDirectory2):
                                if int(os.path.getmtime(winfileOrDirectory2)) != int(member.mtime):
                                    logging.info("'{0}' not up to date because of date change in '{1}'. ({2})".format(tgzFile, winfileOrDirectory2, howToIgnoreString))
                                    return False;
        return True

    @staticmethod
    def WalkWithoutHidden(directory, filter_hidden):
        if not HadoopCopier.is_hidden(directory, filter_hidden):
            for sub in os.listdir(directory):
                subfull = os.path.join(directory,sub)
                if os.path.isfile(subfull):
                    if not HadoopCopier.is_hidden(subfull, filter_hidden):
                        yield subfull
                else:
                    for file in HadoopCopier.WalkWithoutHidden(subfull, filter_hidden):
                        yield file


    @staticmethod
    def dont_filter(tarinfo, directory1):
        filenamein = tarinfo.name.replace("/",os.path.sep)  #"/" is the tar separator
        filename = os.path.join(directory1,filenamein)
        logging.info("Adding file to tar '{0}'".format(filename))
        return tarinfo
    
    @staticmethod
    def tarfile_filter_hidden(tarinfo, directory1, filter_hidden):
        filenamein = tarinfo.name.replace("/",os.path.sep)  #"/" is the tar separator
        filename = os.path.join(directory1,filenamein)
        if HadoopCopier.is_hidden(filename, filter_hidden):
            #logging.info("skipping '{0}'".format(filenamein))
            return None
        else:
            logging.info("Adding file to tar '{0}'".format(filename))
            return tarinfo
        
    @staticmethod
    def is_hidden(filepath, filter_hidden):
        if not filter_hidden:
            return False
        else:
            name = os.path.basename(os.path.abspath(filepath))
            return name.startswith('.') or HadoopCopier.has_hidden_attribute(filepath) or HadoopCopier.containsDotIgnoreTgzChange(filepath)

    _ignoreTgzChangeFileName = ".ignoreTgzChange"

    #!! test that this stops it from look down below
    @staticmethod
    def containsDotIgnoreTgzChange(filepath):
        signalPath = os.path.join(filepath,HadoopCopier._ignoreTgzChangeFileName)
        return os.path.exists(signalPath)

    @staticmethod
    def has_hidden_attribute(filepath):
        try:
            attrs = ctypes.windll.kernel32.GetFileAttributesW(unicode(filepath))
            assert attrs != -1
            result = bool(attrs & 2)
        except (AttributeError, AssertionError):
            result = False
        return result

    def output(self,item):
        outFileList = Hadoop.RecursivelyGetAllOutputs(item)
        for outfile in outFileList:
            if os.path.exists(outfile):
                os.remove(outfile)
            hdfsOutFile = self.remotewd + "/" + outfile
            subprocess.check_output("%HADOOP_HOME%\\bin\\hadoop fs -copyToLocal {0} {1}\n".format(hdfsOutFile, outfile),stderr=subprocess.STDOUT,shell=True)
