
from pysnptools.util.mapreduce1.runner import *
import os
import subprocess, sys, os.path
import multiprocessing
import pysnptools.util as pstutil
import pdb
import logging
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import cPickle as pickle

class HPC(Runner):
    '''
    Old code to run on a Microsoft Widows HPC Cluster. Not currently supported.
    '''
    #!!LATER make it (and Hadoop) work from root directories -- or give a clear error message
    def __init__(self, taskcount, clustername, fileshare, priority="Normal", unit="core", mkl_num_threads=None, runtime="infinite", remote_python_parent=None,
                update_remote_python_parent=False, min=None, max=None, excluded_nodes=[], template=None, nodegroups=None, skipinputcopy=False, node_local=True,clean_up=True,preemptable=True,FailOnTaskFailure=False,logging_handler=logging.StreamHandler(sys.stdout)):
        logger = logging.getLogger()
        if not logger.handlers:
            logger.setLevel(logging.INFO)
        for h in list(logger.handlers):
            logger.removeHandler(h)
        logger.addHandler(logging_handler)
        if logger.level == logging.NOTSET:
            logger.setLevel(logging.INFO)

        self.taskcount = taskcount
        self.clustername = clustername
        self.fileshare = fileshare

        self.priority = priority
        self.runtime = runtime
        self.unit = unit
        self.excluded_nodes = excluded_nodes
        self.min = min
        self.max = max
        self.remote_python_parent = remote_python_parent
        self.update_remote_python_parent = update_remote_python_parent
        self.CheckUnitAndMKLNumThreads(mkl_num_threads, unit)
        self.skipinputcopy=skipinputcopy
        self.template = template
        self.nodegroups = nodegroups
        self.node_local = node_local
        self.clean_up = clean_up
        self.preemptable = preemptable
        self.FailOnTaskFailure = FailOnTaskFailure
      
    def run(self, distributable):
        # Check that the local machine has python path set
        localpythonpath = os.environ.get("PYTHONPATH")#!!should it be able to work without pythonpath being set (e.g. if there was just one file)? Also, is None really the return or is it an exception.
        if localpythonpath is None: raise Exception("Expect local machine to have 'pythonpath' set")

        remotepythoninstall = self.check_remote_pythoninstall()

        remotewd, run_dir_abs, run_dir_rel, nodelocalwd = self.create_run_dir()
        pstutil.create_directory_if_necessary(os.path.join(remotewd, distributable.tempdirectory), isfile=False) #create temp directory now so that cluster tasks won't try to create it many times at once
        result_remote = os.path.join(run_dir_abs,"result.p")

        self.copy_python_settings(run_dir_abs)

        inputOutputCopier = HPCCopier(remotewd,skipinput=self.skipinputcopy) #Create the object that copies input and output files to where they are needed

        inputOutputCopier.input(distributable) # copy of the input files to where they are needed (i.e. the cluster)

        remotepythonpath = self.FindOrCreateRemotePythonPath(localpythonpath, run_dir_abs)

        batfilename_rel = self.create_bat_file(distributable, remotepythoninstall, remotepythonpath, remotewd, run_dir_abs, run_dir_rel, result_remote, nodelocalwd, distributable)

        self.submit_to_cluster(batfilename_rel, distributable, remotewd, run_dir_abs, run_dir_rel, nodelocalwd)

        inputOutputCopier.output(distributable) # copy the output file from where they were created (i.e. the cluster) to the local computer

        assert os.path.exists(result_remote), "The HPC job produced no result (and, thus, likely failed)"
        with open(result_remote, mode='rb') as f:
            result = pickle.load(f)

        #logging.info('Done: HPC runner is running a distributable. Returns {0}'.format(result))
        return result



    def CheckUnitAndMKLNumThreads(self, mkl_num_threads, unit):
        if unit.lower() == "core":
            if mkl_num_threads is not None and mkl_num_threads!=1 : raise Exception("When 'unit' is 'core', mkl_num_threads must be unspecified or 1")
            self.mkl_num_threads = 1
        elif unit.lower() == "socket":
            if mkl_num_threads is None : raise Exception("When 'unit' is 'socket', mkl_num_threads must be specified")
            self.mkl_num_threads = mkl_num_threads
        elif unit.lower() == "node":
            self.mkl_num_threads = mkl_num_threads
        else :
            raise Exception("Expect 'unit' to be 'core', 'socket', or 'node'")

    def copy_python_settings(self, run_dir_abs):
        #localuserprofile = os.environ.get("USERPROFILE")
        user_python_settings=".continuum"
        python_settings=os.path.join(self.fileshare,user_python_settings)
        if os.path.exists(python_settings):
            import shutil
            remote_user_python_settings=os.path.join(run_dir_abs,user_python_settings)
            shutil.copytree(python_settings,remote_user_python_settings)


    def FindOrCreateRemotePythonPath(self, localpythonpath, run_dir_abs):
        if self.remote_python_parent is None:
            remotepythonpath = self.CopySource(localpythonpath, run_dir_abs)
        else:
            pstutil.create_directory_if_necessary(self.remote_python_parent,isfile=False)
            list = []
            for rel in os.listdir(self.remote_python_parent):
                list.append(os.path.join(self.remote_python_parent,rel))
            remotepythonpath = ";".join(list)
            if self.update_remote_python_parent:
                remotepythonpath = self.CopySource(localpythonpath, run_dir_abs)
        
        return remotepythonpath

    def numString(self):
        if self.min is None and self.max is None:
            return " -Num{0} *-*".format(self.unit.capitalize())
        if self.min is None:
            return " -Num{0} {1}".format(self.unit.capitalize(), self.max)
        if self.max is None:
            return " -Num{0} {1}-*".format(self.unit.capitalize(), self.min)
        return " -Num{0} {1}-{2}".format(self.unit.capitalize(), self.min, self.max)

    def submit_to_cluster(self, batfilename_rel, distributable, remotewd, run_dir_abs, run_dir_rel, nodelocalwd):
        stdout_dir_rel = os.path.join(run_dir_rel,"stdout")
        stdout_dir_abs = os.path.join(run_dir_abs,"stdout")
        pstutil.create_directory_if_necessary(stdout_dir_abs, isfile=False)
        stderr_dir_rel = os.path.join(run_dir_rel,"stderr")
        stderr_dir_abs = os.path.join(run_dir_abs,"stderr")
        pstutil.create_directory_if_necessary(stderr_dir_abs, isfile=False)
        
        if len(self.excluded_nodes) > 0:
            excluded_nodes = "Set-HpcJob -Id $r.Id -addExcludedNodes {0}".format(", ".join(self.excluded_nodes))
        else:
            excluded_nodes = ""


        #create the Powershell file
        psfilename_rel = os.path.join(run_dir_rel,"dist.ps1")
        psfilename_abs = os.path.join(run_dir_abs,"dist.ps1")
        pstutil.create_directory_if_necessary(psfilename_abs, isfile=True)
        with open(psfilename_abs, "w") as psfile:
            psfile.write(r"""Add-PsSnapin Microsoft.HPC
        Set-Content Env:CCP_SCHEDULER {0}
        $r = New-HpcJob -Name "{7}" -Priority {8}{12}{14}{16} -RunTime {15} -FailOnTaskFailure {23} #-Preemptable {22}
        $r.Id
        if ({20})
        {10}
            $from = "{4}"
            $to = "{17}"
            Add-HpcTask -Name NodePrep    -JobId $r.Id -Type NodePrep                -CommandLine "${{from}}\{18}"        -StdOut "${{from}}\{2}\nodeprep.txt"    -StdErr "${{from}}\{3}\nodeprep.txt"    -WorkDir .
            Add-HpcTask -Name Parametric  -JobId $r.Id -Parametric -Start 0 -End {1} -CommandLine "${{from}}\{6} * {5}"   -StdOut "${{from}}\{2}\*.txt"    -StdErr "${{from}}\{3}\*.txt"                  -WorkDir $to
            Add-HpcTask -Name Reduce      -JobId $r.Id -Depend Parametric            -CommandLine "${{from}}\{6} {5} {5}" -StdOut "${{from}}\{2}\reduce.txt"      -StdErr "${{from}}\{3}\reduce.txt"      -WorkDir $to
            {21}Add-HpcTask -Name NodeRelease -JobId $r.Id -Type NodeRelease         -CommandLine "${{from}}\{19}"        -StdOut "${{from}}\{2}\noderelease.txt" -StdErr "${{from}}\{3}\noderelease.txt" -WorkDir .
        {11}
        else
        {10}
            Add-HpcTask -Name Parametric -JobId $r.Id -Parametric -Start 0 -End {1} -CommandLine "{6} * {5}" -StdOut "{2}\*.txt" -StdErr "{3}\*.txt" -WorkDir {4}
            Add-HpcTask -Name Reduce -JobId $r.Id -Depend Parametric -CommandLine "{6} {5} {5}" -StdOut "{2}\reduce.txt" -StdErr "{3}\reduce.txt" -WorkDir {4}
        {11}

        {13}
        Submit-HpcJob -Id $r.Id
        $j = Get-HpcJob -Id $r.Id
        $i = $r.id
        $s = 10

        while(($j.State -ne "Finished") -and ($j.State -ne "Failed") -and ($j.State -ne "Canceled"))
        {10}
            $x = $j.State
            Write-Host "${10}x{11}. Job# ${10}i{11} sleeping for ${10}s{11}"
            Start-Sleep -s $s
            if ($s -ge 60)
            {10}
            $s = 60
            {11}
            else
            {10}
                $s = $s * 1.1
            {11}
           $j.Refresh()
        {11}

        """                 .format(
                                self.clustername,   #0
                                self.taskcount-1,   #1
                                stdout_dir_rel,     #2
                                stderr_dir_rel,     #3
                                remotewd,           #4 fileshare wd
                                self.taskcount,     #5
                                batfilename_rel,    #6
                                self.maxlen(str(distributable),50),      #7
                                self.priority,      #8
                                self.unit,          #9 -- not used anymore,. Instead #12 sets unit
                                "{",                #10
                                "}",                #11
                                self.numString(),   #12
                                excluded_nodes,     #13
                                ' -templateName "{0}"'.format(self.template) if self.template is not None else "", #14
                                self.runtime,       #15 RuntimeSeconds
                                ' -NodeGroups "{0}"'.format(self.nodegroups) if self.nodegroups is not None else "", #16
                                nodelocalwd,        #17 the node-local wd
                                batfilename_rel[0:-8]+"nodeprep.bat", #18
                                batfilename_rel[0:-8]+"noderelease.bat", #19
                                1 if self.node_local else 0,             #20
                                "",                                      #21 always run release task
                                self.preemptable,                        #22
                                '$true' if self.FailOnTaskFailure else '$false',   #23
                                ))
        assert batfilename_rel[-8:] == "dist.bat", "real assert"
        import subprocess
        proc = subprocess.Popen(["powershell.exe", "-ExecutionPolicy", "Unrestricted", psfilename_abs], cwd=os.getcwd())
        if not 0 == proc.wait(): raise Exception("Running powershell cluster submit script results in non-zero return code")

    #move to utils?
    @staticmethod
    def maxlen(s,max):
        '''
        Truncate cluster job name if longer than max.
        '''
        if len(s) <= max:
            return s
        else:
            #return s[0:max-1]
            return s[-max:]  #JL: I prefer the end of the name rather than the start

    
    def create_distributablep(self, distributable, run_dir_abs, run_dir_rel):
        distributablep_filename_rel = os.path.join(run_dir_rel, "distributable.p")
        distributablep_filename_abs = os.path.join(run_dir_abs, "distributable.p")
        with open(distributablep_filename_abs, mode='wb') as f:
            pickle.dump(distributable, f, pickle.HIGHEST_PROTOCOL)
        return distributablep_filename_rel, distributablep_filename_abs

    @staticmethod
    def FindDirectoriesToExclude(localpythonpathdir):
        logging.info("Looking in '{0}' for directories to skip".format(localpythonpathdir))
        xd_string = " /XD $TF /XD .git"
        for root, dir, files in os.walk(localpythonpathdir):
            for file in files:
             if file.lower() == ".ignoretgzchange":
                 xd_string += " /XD {0}".format(root)
        return xd_string

    def CopySource(self,localpythonpath, run_dir_abs):
        
        if self.update_remote_python_parent:
            remote_python_parent = self.remote_python_parent
        else:
            remote_python_parent = run_dir_abs + os.path.sep + "pythonpath"
        pstutil.create_directory_if_necessary(remote_python_parent, isfile=False)
        remotepythonpath_list = []
        for i, localpythonpathdir in enumerate(localpythonpath.split(';')):
            remotepythonpathdir = os.path.join(remote_python_parent, str(i))
            remotepythonpath_list.append(remotepythonpathdir)
            xd_string = HPC.FindDirectoriesToExclude(localpythonpathdir)
            xcopycommand = 'robocopy /s {0} {1}{2}'.format(localpythonpathdir,remotepythonpathdir,xd_string)
            logging.info(xcopycommand)
            os.system(xcopycommand)

        remotepythonpath = ";".join(remotepythonpath_list)
        return remotepythonpath

    def create_bat_file(self, distributable, remotepythoninstall, remotepythonpath, remotewd, run_dir_abs, run_dir_rel, result_remote, nodelocalwd, create_bat_file):
        path_share_list = [r"",r"Scripts"]
        remotepath_list = []
        for path_share in path_share_list:
            path_share_abs = os.path.join(remotepythoninstall,path_share)
            if not os.path.isdir(path_share_abs): raise Exception("Expect path directory at '{0}'".format(path_share_abs))
            remotepath_list.append(path_share_abs)
        remotepath = ";".join(remotepath_list)

        distributablep_filename_rel, distributablep_filename_abs = self.create_distributablep(distributable, run_dir_abs, run_dir_rel)

        distributable_py_file = os.path.join(os.path.dirname(__file__),"..","distributable.py")
        if not os.path.exists(distributable_py_file): raise Exception("Expect file at " + distributable_py_file + ", but it doesn't exist.")
        localfilepath, file = os.path.split(distributable_py_file)

        for remote_path_part in remotepythonpath.split(';'):
            remoteexe = os.path.join(remote_path_part,"fastlmm","util",file)
            if os.path.exists(remoteexe):
                break #not continue
            remoteexe = None
        assert remoteexe is not None, "Could not find '{0}' on remote python path. Is fastlmm on your local python path?".format(file)

        #run_dir_rel + os.path.sep + "pythonpath" + os.path.sep + os.path.splitdrive(localfilepath)[1]

        #result_remote2 = result_remote.encode("string-escape")
        command_string = remoteexe + r""" "{0}" """.format(distributablep_filename_abs) + r""" "LocalInParts(%1,{0},mkl_num_threads={1},result_file=""{2}"",run_dir=""{3}"") " """.format(
            self.taskcount,
            self.mkl_num_threads,
            "result.p",
            run_dir_abs.encode("string-escape"))
        batfilename_rel = os.path.join(run_dir_rel,"dist.bat")
        batfilename_abs = os.path.join(run_dir_abs,"dist.bat")
        pstutil.create_directory_if_necessary(batfilename_abs, isfile=True)
        matplotlibfilename_rel = os.path.join(run_dir_rel,".matplotlib")
        matplotlibfilename_abs = os.path.join(run_dir_abs,".matplotlib")
        pstutil.create_directory_if_necessary(matplotlibfilename_abs, isfile=False)
        pstutil.create_directory_if_necessary(matplotlibfilename_abs + "/tex.cache", isfile=False)
        ipythondir_rel = os.path.join(run_dir_rel,".ipython")
        ipythondir_abs = os.path.join(run_dir_abs,".ipython")
        pstutil.create_directory_if_necessary(ipythondir_abs, isfile=False)
        with open(batfilename_abs, "w") as batfile:
            batfile.write("set path={0};%path%\n".format(remotepath))
            batfile.write("set PYTHONPATH={0}\n".format(remotepythonpath))
            batfile.write("set USERPROFILE={0}\n".format(run_dir_abs))
            batfile.write("set MPLCONFIGDIR={0}\n".format(matplotlibfilename_abs))
            batfile.write("set IPYTHONDIR={0}\n".format(ipythondir_abs))
            batfile.write("python {0}\n".format(command_string))

        if (self.node_local):
            with open( os.path.join(run_dir_abs,"nodeprep.bat"), "w") as prepfile:
                prepfile.write(r"""set f="{0}"{1}""".format(remotewd,'\n'))
                prepfile.write(r"""set t="{0}"{1}""".format(nodelocalwd,'\n'))
                prepfile.write("if not exist %t% mkdir %t%\n")
                with open( os.path.join(run_dir_abs,"noderelease.bat"), "w") as releasefile:
                    releasefile.write(r"""set f="{0}"{1}""".format(remotewd,'\n'))
                    releasefile.write(r"""set t="{0}"{1}""".format(nodelocalwd,'\n'))
                    inputOutputCopier = HPCCopierNodeLocal(prepfile,releasefile,self.clean_up) #Create the object that copies input and output files to where they are needed
                    inputOutputCopier.input(distributable) # copy of the input files to where they are needed (i.e. to the cluster)
                    inputOutputCopier.output(distributable) # copy of the output files to where they are needed (i.e. off the cluster)
                    releasefile.write("rmdir /s %t%\n")
                    releasefile.write("exit /b 0\n")


        return batfilename_rel

    def check_remote_pythoninstall(self):
        remotepythoninstall = r"\\GCR\Scratch\RR1\escience\pythonInstallD"  #!!! don't hardwire this
        if not os.path.isdir(remotepythoninstall): raise Exception("Expect Python and related directories at '{0}'".format(remotepythoninstall))

        return remotepythoninstall

    def create_run_dir(self):
        username = os.environ["USERNAME"]
        localwd = os.getcwd()
        #!!make an option to specify the full remote WD. Also what is the "\\\\" case for?
        if localwd.startswith("\\\\"):
            remotewd = self.fileshare + os.path.sep + username +os.path.sep + "\\".join(localwd.split('\\')[4:])
            nodelocalwd =  "d:\scratch\escience" + os.path.sep + username +os.path.sep + "\\".join(localwd.split('\\')[4:]) #!!!const
        else:
            remotewd = self.fileshare + os.path.sep + username + os.path.splitdrive(localwd)[1]  #using '+' because 'os.path.join' isn't work with shares
            nodelocalwd = "d:\scratch\escience" + os.path.sep + username + os.path.splitdrive(localwd)[1]  #!!! const
        import datetime
        now = datetime.datetime.now()
        run_dir_rel = os.path.join("runs",pstutil._datestamp(appendrandom=True))
        run_dir_abs = os.path.join(remotewd,run_dir_rel)
        pstutil.create_directory_if_necessary(run_dir_abs,isfile=False)


        return remotewd, run_dir_abs, run_dir_rel, nodelocalwd


class HPCCopier(object): #Implements ICopier

    def __init__(self, remotewd, skipinput=False):
        self.remotewd = remotewd
        self.skipinput=skipinput

    def input(self,item):
        if self.skipinput:
            return
        if isinstance(item, str):
            itemnorm = os.path.normpath(item)
            remote_file_name = os.path.join(self.remotewd,itemnorm)
            remote_dir_name,ignore = os.path.split(remote_file_name)
            pstutil.create_directory_if_necessary(remote_file_name)
            xcopycommand = "xcopy /d /e /s /c /h /y {0} {1}".format(itemnorm, remote_dir_name)
            logging.info(xcopycommand)
            rc = os.system(xcopycommand)
            print "rc=" +str(rc)
            if rc!=0: raise Exception("xcopy cmd failed with return value={0}, from cmd {1}".format(rc,xcopycommand))
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        # else -- do nothing

    def output(self,item):
        if isinstance(item, str):
            itemnorm = os.path.normpath(item)
            pstutil.create_directory_if_necessary(itemnorm)
            remote_file_name = os.path.join(self.remotewd,itemnorm)
            local_dir_name,ignore = os.path.split(itemnorm)
            assert os.path.exists(remote_file_name), "Don't see expected file '{0}'. Did the HPC job fail?".format(remote_file_name)
            #xcopycommand = "xcopy /d /e /s /c /h /y {0} {1}".format(remote_file_name, local_dir_name) # we copy to the local dir instead of the local file so that xcopy won't ask 'file or dir?'
            xcopycommand = "xcopy /d /c /y {0} {1}".format(remote_file_name, local_dir_name) # we copy to the local 
            logging.info(xcopycommand)
            rc = os.system(xcopycommand)
            if rc!=0: logging.info("xcopy cmd failed with return value={0}, from cmd {1}".format(rc,xcopycommand))
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        # else -- do nothing

class HPCCopierNodeLocal(object): #Implements ICopier

    def __init__(self, fileprep, filerelease, clean_up):
        self.fileprep = fileprep
        self.filerelease = filerelease
        self.clean_up = clean_up

    def input(self,item):
        if isinstance(item, str):
            itemnorm = os.path.normpath(item)
            dirname = os.path.dirname(itemnorm)
            self.fileprep.write("if not exist %t%\{0} mkdir %t%\{0}\n".format(dirname))
            self.fileprep.write("xcopy /d /e /s /c /h /y %f%\{0} %t%\{1}\\\n".format(itemnorm,dirname))
            if self.clean_up:
                self.filerelease.write("del %t%\{0}\n".format(itemnorm))
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        # else -- do nothing

    def output(self,item):
        if isinstance(item, str):
            itemnorm = os.path.normpath(item)
            dirname = os.path.dirname(itemnorm)
            self.filerelease.write("xcopy /d /e /s /c /h /y %t%\{0} %f%\{1}\\\n".format(itemnorm,dirname))
            if self.clean_up:
                self.filerelease.write("del %t%\{0}\n".format(itemnorm))
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        # else -- do nothing
