"""Tools for reading and writing files, locally or across clusters.
"""

import os
import threading
import re
import time
import logging
import random


def ip_address():
    '''
    Return the ip address of this computer.
    '''
    import socket
    #see http://stackoverflow.com/questions/166506/finding-local-ip-addresses-using-pythons-stdlib
    return ([l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1],[[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]) 

def ip_address_pid():
    '''
    Return the ip address of this computer and process id of the executing process.
    '''
    ip = ip_address()
    pid = os.getpid()
    tid = threading.current_thread().ident
    return "{0}.{1}.{2}".format(ip,pid,tid)

class _DibLib(object):
    def __init__(self,unique_name,dib_path,dir_path,dib_name):
        self.dib_path=dib_path
        self.dir_path=dir_path
        self.dibs_fn = "{0}_{1}.txt".format(dib_name,unique_name)
        self.dibs_pattern = re.compile("^{0}_.*\.txt$".format(dib_name))
        if dib_name == "dibs": #could make two subclasses, but this seems simpler
            self.test = self.tree_copy
        elif dib_name == "azure_storage_dib":
            self.test = self.only_first_goes
        else:
            raise Exception("Don't know what test goes with dib_name of '{0}'".format(dib_name))

    def tree_copy(self,priority,dibs_time,try_index):
        copy_count = len([f for f in self.dir_path.walk() if PeerToPeer._copy_main_pattern.match(f)]) #count main + copies
        logging.info("PeerToPeer priority is {0} and copy_count is {1}".format(priority,copy_count))
        if priority+1 <= copy_count * 2:
            return 'go'
        else:
            return 'wait'
    
    def only_first_goes(self,priority,dibs_time,try_index):
        logging.info("P2P repair priority is {0}".format(priority))
        if priority == 0 and try_index == 0:
            return "azure" #If dibs first, then go
        if priority == 0:
            logging.info("Priority was > 0, but there is now a p2p file to copy, so the wait is over")
            return "fixed"
        for f in self.dir_path.walk():
            if PeerToPeer._copy_main_pattern.match(f) and self.dir_path._simple_getmtime(f) > dibs_time:
                logging.info("Priority > 0, but there is now a p2p file to copy, so the wait is over")
                return "fixed" #If dibs later, but there is a new "main" or "copy" file, then go
        else:
            return "wait"  # else, wait

    def ever(self):
        i = 0
        while True:
            yield i
            i += 1

    def wait_for_turn(self):
        self.dib_path.save(self.dibs_fn,"")
        time.sleep(1.0) # Sleep one second to (hopefully) see dibs made at exactly the same time.
        dibs_time = self.dib_path.getmtime(self.dibs_fn)
        sleep_time = 5.0
        for try_index in self.ever():
            priority = self.get_priority(dibs_time) # Count number of remaining dibs files that were created before mine. That is my priority.
            status = self.test(priority,dibs_time,try_index)
            if status != 'wait':
                logging.info("I can go! (status='{0}')".format(status))
                return status
            logging.info("sleep for {0}".format(sleep_time))
            time.sleep(sleep_time)
            sleep_time = min(sleep_time*random.uniform(1.0,1.2),60) #wait 0 to 10% longer, but no more than a minute

    def remove_dibs(self):
        if self.dib_path._simple_file_exists(self.dibs_fn):
            self.dib_path._simple_remove(self.dibs_fn)

    def get_priority(self,dibs_time):
        priority = 0
        sequence = (file for file in self.dib_path.walk() if self.dibs_pattern.match(file))
        for other_dib in sequence:
            try: #It might disappear out from under us
                other_time = self.dib_path.getmtime(other_dib)
            except:
                other_time = None
            if  other_time is None or other_time < dibs_time or (other_time==dibs_time and other_dib < self.dibs_fn):
                priority += 1
        return priority

 

from pysnptools.util.filecache.filecache import FileCache
from pysnptools.util.filecache.localcache import LocalCache
from pysnptools.util.filecache.peertopeer import PeerToPeer


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
